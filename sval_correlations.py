import os
import copy
import pdb
import pickle

import numpy as np
import matplotlib.pyplot as plt
import corner
import radvel
import radvel.fitting
from radvel.plot import orbit_plots

"""Test code to measure correlations between each individual Keplerian
model and S-values.

"""

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

        # Separate the RVs into independent instrument datasets.
        self.rv_dict    = {}
        self.times_dict = {}

        if data is not None:
            self.times = data['jd']
            for tel in self.tels:
                self.rv_dict[tel]    = data.query('tel == "{}"'.format(tel)).mnvel
                self.times_dict[tel] = data.query('tel == "{}"'.format(tel)).jd
        else:
            self.times = self.post.likelihood.x
            for like in self.post.likelihood.like_list:
                tel = like.telvec[0]
                self.rv_dict[tel]    = like.y
                self.times_dict[tel] = like.x

        self.rv_resid_dict = {}
        self.starname = starname
        self.outdir   = outdir
        full_outdir = os.path.join(os.getcwd(), outdir)
        if not os.path.exists(full_outdir):
            os.mkdir(full_outdir)

    def residuals(self, tel, num_model):
        """Compute residuals for one Keplerian model, in one dataset.

        """
        rvs = copy.deepcopy(self.rv_dict[tel])
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

    def correlate(self, tel, num_model):
        """Measure correlation for one Keplerian model, in one dataset.

        """
        pass

    def make_all_residuals(self):
        for tel in self.tels:
            for n in np.arange(1, self.nplanets+1):
                self.residuals(tel, n)

    def stack_one_tel(self, tel):
        """Make stacked correlation plot for one dataset.

        """
        sval = self.svals[tel]
        num_plots = self.nplanets
        figwidth=7.5
        fig = plt.figure(figsize=(figwidth, 0.5*figwidth*(num_plots+1)))
        plt.title('{}, {}'.format(self.starname, tel))

        for i in np.arange(1,num_plots+1):
            ax = fig.add_subplot(num_plots, 1, i)
            ax.set_xlim([np.amin(sval)-0.005, np.amax(sval)+0.005])
            ax.scatter(sval, self.rv_resid_dict[tel+str(i)], c='black')
            ax.set_ylabel(r'$\Delta$RV$_{}$'.format(i))
            #ax.legend(loc=0)
            if i != num_plots:
                ax.tick_params(axis='x', which='both', direction='in',
                               bottom='on', top='off', labelbottom='off')
            plt.subplots_adjust(hspace=0.00)
            plt.subplots_adjust(wspace=0.00)

        ax.set_xlabel('S-value')

        fig.savefig(self.outdir+'/{}_{}.pdf'.format(self.starname, tel))
