import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils.io

import rvsearch
from rvsearch import utils

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'HIP52942A'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.05, mstar_err=.06)

# load in data
data = cpsutils.io.loadcps('hip52942a', hires_rk=True, hires_rj=True,
                           ctslim=3000, binsize=0.5)
#data = data[data['obnm'] != 'rj81.412']  # low counts, only 30k and also poor seeing
data['time'] = data['bjd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=30477.0)
    params['tc1'] = radvel.Parameter(value=2451703.3)
    params['k1'] = radvel.Parameter(value=4252.0)
    params['secosw1'] = radvel.Parameter(value=0.23)
    params['sesinw1'] = radvel.Parameter(value=0.65)
    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-20.0, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=7.0, vary=True)
params['gamma_k'] = radvel.Parameter(value=0.0, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=10.0, vary=True)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 30.0),
    radvel.prior.HardBounds('jit_j', 0.0, 30.0),
]
