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

starname = 'HD 239960'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
# stellar = dict(mstar=0.57, mstar_err=.02)

# load in data
data = cpsutils.io.loadcps('239960', hires_rk=True, hires_rj=True, lick=False, ctslim=3000, binsize=0.5)
data = data[data['obnm'] != 'rj179.332']  # low counts, only 30k and also poor seeing

data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=19099.6)
    params['tc1'] = radvel.Parameter(value=2453800.0)
    params['k1'] = radvel.Parameter(value=648.7)
    params['e1'] = radvel.Parameter(value=0.44)
    params['w1'] = radvel.Parameter(value=-2.74)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-3963.8, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=9.6, vary=True)
params['gamma_k'] = radvel.Parameter(value=15698.8, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=4.7, vary=True)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.1, 50.0),
    radvel.prior.HardBounds('jit_j', 0.1, 50.0),
]
