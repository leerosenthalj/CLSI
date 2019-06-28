import pdb

import pandas as pd
import numpy as np
import radvel

import rvsearch
from rvsearch import utils

starname = 'HD125612'
nplanets = 3
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('./setup_data/vst125612.csv')
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=556.29)
    params['tp1'] = radvel.Parameter(value=2463233.1)
    params['k1'] = radvel.Parameter(value=80.74)
    params['e1'] = radvel.Parameter(value=0.48)
    params['w1'] = radvel.Parameter(value=0.72)

    params['per2'] = radvel.Parameter(value=4.15458)
    params['tp2'] = radvel.Parameter(value=2463058.4)
    params['k2'] = radvel.Parameter(value=6.0)
    params['e2'] = radvel.Parameter(value=0.21)
    params['w2'] = radvel.Parameter(value=-2.01)

    params['per3'] = radvel.Parameter(value=2859.01)
    params['tp3'] = radvel.Parameter(value=2463138.98)
    params['k3'] = radvel.Parameter(value=100.16)
    params['e3'] = radvel.Parameter(value=0.13)
    params['w3'] = radvel.Parameter(value=-0.68)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=0)
params['jit_j'] = radvel.Parameter(value=1.0)
params['gamma_k'] = radvel.Parameter(value=0)
params['jit_k'] = radvel.Parameter(value=1.0)
#params['gamma_lick'] = radvel.Parameter(value=222.728)
#params['jit_lick'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]