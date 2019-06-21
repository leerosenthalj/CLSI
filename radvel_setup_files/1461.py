import pdb

import pandas as pd
import numpy as np
import radvel

import rvsearch
from rvsearch import utils

starname = 'HD1461'
nkeplers = 3
instnames = ['k', 'j', 'apf']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('vst1461.csv')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nkeplers, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=5.77132)
    params['tp1'] = radvel.Parameter(value=7.14416)
    params['k1'] = radvel.Parameter(value=2.57138)
    params['e1'] = radvel.Parameter(value=0.093)
    params['w1'] = radvel.Parameter(value=-0.436)

    params['per2'] = radvel.Parameter(value=13.5053)
    params['tp2'] = radvel.Parameter(value=5205.13)
    params['k2'] = radvel.Parameter(value=2.03808)
    params['e2'] = radvel.Parameter(value=0.175)
    params['w2'] = radvel.Parameter(value=-2.249)

    params['per3'] = radvel.Parameter(value=3941.96)
    params['tp3'] = radvel.Parameter(value=3796.67)
    params['k3'] = radvel.Parameter(value=2.10581)
    params['e3'] = radvel.Parameter(value=0.080)
    params['w3'] = radvel.Parameter(value=0.360)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=0)
params['jit_j'] = radvel.Parameter(value=2.)
params['gamma_k'] = radvel.Parameter(value=0)
params['jit_k'] = radvel.Parameter(value=2.)
params['gamma_apf'] = radvel.Parameter(value=0)
params['jit_apf'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nkeplers ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0)
]
