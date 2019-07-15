import pdb

import pandas as pd
import numpy as np
import radvel
import cpsutils
import cpsutils.io

import rvsearch
from rvsearch import utils

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'HD 156279'
nplanets = 2
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=0.96, mstar_err=.05)

# load in data
data = cpsutils.io.loadcps('156279', hires_rk=True, hires_rj=True,
                           lick=False, ctslim=3000, binsize=2.0)
data['tel'] = data['tel'].str.decode('utf-8')
data['time'] = data['jd']
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nplanets,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=133.39)
    params['tc1'] = radvel.Parameter(value=2455464.56)
    params['k1'] = radvel.Parameter(value=509.7)
    params['e1'] = radvel.Parameter(value=0.64)
    params['w1'] = radvel.Parameter(value=-1.64)
    params['per2'] = radvel.Parameter(value=4722.1)
    params['tc2'] = radvel.Parameter(value=2455838.4)
    params['k2'] = radvel.Parameter(value=121.8)
    params['e2'] = radvel.Parameter(value=0.25)
    params['w2'] = radvel.Parameter(value=1.77)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-60, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=2.15, vary=True)
params['gamma_k'] = radvel.Parameter(value=-60, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=2.72, vary=True)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.1, 50.0),
    radvel.prior.HardBounds('jit_j', 0.1, 50.0),
]
