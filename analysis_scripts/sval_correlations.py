import os
import copy
import pdb
import pickle
from math import floor, log10

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec
from matplotlib.ticker import ScalarFormatter, NullFormatter

import corner
import radvel
import radvel.fitting
from radvel.plot import orbit_plots

"""Test code to measure correlations between each individual Keplerian
model and S-values.

"""

# Define function to round to significant figures.
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1) #log10(abs(x))
round_err_vec = np.vectorize(round_sig)
def round_to_ref(x, y):
    return round(x, -int(floor(log10(y))))
round_vec = np.vectorize(round_to_ref)

class Correlator(object):
    """Class to measure correlations between each individual Keplerian
    model and S-values.

    Args:
        post (radvel.Posterior): Optional posterior with known planet params.
        svals (dict): Dictionary of S-value arrays for separate instruments.
            NOTE: svals instrument keys must match posterior keys.
        starname (str): String, used to name the output directory.

    """

    def __init__(self, post, svals, data=None, starname='star', outdir='correlate'):

        self.post     = post
        self.nplanets = post.params.num_planets

        self.svals = svals
        self.tels  = svals.keys() # NOT including instruments without Svals.

        self.linfit = {}

        # Separate the RVs into independent instrument datasets.
        self.mnvel_dict  = {}
        self.errvel_dict = {}
        self.times_dict  = {}

        if data is not None:
            self.times = data['jd']
            for tel in self.tels:
                self.mnvel_dict[tel]  = data.query('tel == "{}"'.format(tel)).mnvel
                self.errvel_dict[tel] = data.query('tel == "{}"'.format(tel)).errvel
                self.times_dict[tel]  = data.query('tel == "{}"'.format(tel)).jd
        else:
            self.times = self.post.likelihood.x
            for like in self.post.likelihood.like_list:
                tel = like.telvec[0]
                self.mnvel_dict[tel] = like.y
                self.times_dict[tel] = like.x

        self.rv_resid_dict = {}
        self.starname = starname
        self.outdir   = outdir
        full_outdir   = os.path.join(os.getcwd(), outdir)
        if not os.path.exists(full_outdir):
            os.mkdir(full_outdir)

    def residuals(self, tel, num_model):
        """Compute residuals for one Keplerian model, in one dataset.

        """
        rvs = copy.deepcopy(self.mnvel_dict[tel])
        # Subtract offset and trend/curvature terms.
        rvs -= self.post.params['gamma_{}'.format(tel)].value

        if self.post.params['dvdt'].vary == True:
            rvs -= self.post.params['dvdt'].value * (self.times_dict[tel] -
                   self.post.likelihood.model.time_base)

        if self.post.params['curv'].vary == True:
            rvs -= self.post.params['curv'].value * (self.times_dict[tel] -
                   self.post.likelihood.model.time_base)**2

        # Subtract residuals of all Keplerian models except the selected one.
        if self.nplanets > 1:
            indices = np.arange(1, self.nplanets+1)
            for i in indices[indices != num_model]:
                orbel = [self.post.params['per{}'.format(i)].value,
                         self.post.params['tp{}'.format(i)].value,
                         self.post.params['e{}'.format(i)].value,
                         self.post.params['w{}'.format(i)].value,
                         self.post.params['k{}'.format(i)].value]
                mod  = radvel.kepler.rv_drive(np.array(self.times_dict[tel]), orbel)
                rvs -= mod
                # Save residuals to the dictionary.
                self.rv_resid_dict[tel+str(num_model)] = rvs

    def make_all_residuals(self):
        for tel in self.tels:
            for n in np.arange(1, self.nplanets+1):
                self.residuals(tel, n)

    def correlate(self, tel, num_model):
        """Measure correlation for one Keplerian model, in one dataset.

        """
        sval  = self.svals[tel]
        rvs_n = self.rv_resid_dict[tel+str(num_model)]
        errs  = np.sqrt(self.errvel_dict[tel]**2 +
                        self.post.params['jit_{}'.format(tel)].value**2)

        A    = np.vstack((np.ones_like(sval), sval)).T
        C    = np.diag(errs * errs)
        cov  = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
        b, m = np.dot(cov, np.dot(A.T, np.linalg.solve(C, rvs_n)))
        sigb = np.sqrt(cov[0][0])
        sigm = np.sqrt(cov[1][1])

        self.linfit[tel+str(num_model)] = [b, m, sigb, sigm]


    def stack_one_tel(self, tel):
        """Make stacked correlation plot for one dataset.

        """
        sval = self.svals[tel]
        nmodels = self.nplanets
        figwidth=7.5
        fig = plt.figure(figsize=(figwidth, 0.5*figwidth*(nmodels+1)))
        fig.suptitle('{}, {}'.format(self.starname, tel), fontsize=30)

        # Try Gridspec.
        self.ax_list = []
        gridthing = gridspec.GridSpec(nmodels, 1)
        for i in np.arange(nmodels):#np.arange(1,nmodels+1):
            rvs_i  = self.rv_resid_dict[tel+str(i+1)]
            ax = plt.subplot(gridthing[i, 0])
            #if i == 0:
            #    ax.set_title('{}, {}'.format(self.starname, tel), fontsize=30)
            self.ax_list += [ax]
            ax.set_xlim([np.amin(sval[~np.isnan(sval)]) - 0.001,
                         np.amax(sval[~np.isnan(sval)]) + 0.001])
            if tel == 'a' or tel == 'apf':
                color='green'
                mark='d'
            else:
                color='black'
                mark='o'
            ax.scatter(sval, rvs_i, c=color, marker=mark, s=100, alpha=0.5)
            ax.set_ylabel(r'$\Delta$RV$_{}$'.format(i+1), fontsize=20)

            # If correlation has been measured, plot it.
            if tel+str(i+1) in self.linfit.keys():
                m_err = round_sig(self.linfit[tel+str(i+1)][3])
                m_round   = round_to_ref(self.linfit[tel+str(i+1)][1], m_err)
                ax.plot(sval, self.linfit[tel+str(i+1)][0] +
                              self.linfit[tel+str(i+1)][1]*sval,
                              label=r'$m={0} \pm {1}$'.format(m_round, m_err))
                ax.legend()
            if i != nmodels-1:
                ax.tick_params(axis='x', which='both', direction='in',
                               bottom='on', top='off', labelbottom='off')

            plt.subplots_adjust(hspace=0.00)
            plt.subplots_adjust(wspace=0.00)

        ax.set_xlabel('S-value', fontsize=20)

        fig.savefig(self.outdir+'/{}_{}.pdf'.format(self.starname, tel))
