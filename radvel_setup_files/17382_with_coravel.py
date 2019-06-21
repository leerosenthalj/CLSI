import pdb

import pandas as pd
import numpy as np
import radvel
from cpsutils import io

import rvsearch
from rvsearch import utils

# radvel fit -s 120066.py -d 2019-2-19/master
# radvel mcmc -s 120066.py -d 2019-2-19/master --maxGR 1.001 --minsteps 1000 --nsteps 10000 --minpercent 100

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'HD17382'
nplanets = 1
instnames = ['j', 'CORAVEL']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=0.726, mstar_err=0.03)

# load in data
#data = io.loadcps('190406', apf=True, hires_rj=True, hires_rk=True,
#                  lick=True, verbose=False, ctslim=3000, detrend=False, binsize=2.0)
data_folder = './'
data = utils.read_from_csv(data_folder+'17382_with_coravel.csv', binsize=0.5)
data['time'] = data['jd']
#data['tel'] = data['tel'].str.decode('utf-8')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(1, basis='per tc e w k')
    #From Halbwachs et al. 2018
    params['per1'] = radvel.Parameter(value=5557.)
    params['tc1'] = radvel.Parameter(value=2453545.3)
    params['k1'] = radvel.Parameter(value=2822.8)
    params['e1'] = radvel.Parameter(value=0.65)
    params['w1'] = radvel.Parameter(value=1.989)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params


# initialize the orbit parameters and the orbit model
params = initialize_params()
#params['gamma_j'] = radvel.Parameter(value=-50.279, linear=True, vary=False)
params['gamma_j'] = radvel.Parameter(value=-615.8, linear=True, vary=False)
params['jit_j'] = radvel.Parameter(value=13.6)
#params['gamma_CORAVEL'] = radvel.Parameter(value=-50.279, linear=True, vary=False)
params['gamma_CORAVEL'] = radvel.Parameter(value=357.3, linear=True, vary=False)
params['jit_CORAVEL'] = radvel.Parameter(value=0.01, vary=False)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 50.0),
#    radvel.prior.HardBounds('jit_CORAVEL', 0.0, 50.0),
]
