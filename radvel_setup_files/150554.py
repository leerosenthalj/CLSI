import pdb

import pandas as pd
import numpy as np
import radvel

import rvsearch
from rvsearch import utils

starname = 'HD150554'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('./setup_data/vst150554.csv')
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=3690.25)
    params['tp1'] = radvel.Parameter(value=2456336.25)
    params['k1'] = radvel.Parameter(value=980.47)
    params['e1'] = radvel.Parameter(value=0.57)
    params['w1'] = radvel.Parameter(value=0.44)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=17.9)
params['jit_j'] = radvel.Parameter(value=1.5)
params['gamma_k'] = radvel.Parameter(value=27.7)
params['jit_k'] = radvel.Parameter(value=1.5)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]
