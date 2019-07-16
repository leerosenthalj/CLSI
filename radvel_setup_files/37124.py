import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD37124'
nplanets = 3
instnames = ['k', 'j', 'apf']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('./setup_data/vst37124.csv')
#data = cpsutils.io.loadcps('37124', hires_rk=True, hires_rj=True,
#                           apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=154.25)
    params['tp1'] = radvel.Parameter(value=2454935.77)
    params['k1'] = radvel.Parameter(value=28.76)
    params['e1'] = radvel.Parameter(value=0.05415)
    params['w1'] = radvel.Parameter(value=2.28)

    params['per2'] = radvel.Parameter(value=885.265)
    params['tp2'] = radvel.Parameter(value=2454939.34)
    params['k2'] = radvel.Parameter(value=15.3)
    params['e2'] = radvel.Parameter(value=0.12)
    params['w2'] = radvel.Parameter(value=0.76)

    params['per3'] = radvel.Parameter(value=1761.56)
    params['tp3'] = radvel.Parameter(value=2454422.54)
    params['k3'] = radvel.Parameter(value=13.04)
    params['e3'] = radvel.Parameter(value=0.02382)
    params['w3'] = radvel.Parameter(value=0.8725)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=3.904)
params['jit_j'] = radvel.Parameter(value=3.)
params['gamma_k'] = radvel.Parameter(value=-1.107)
params['jit_k'] = radvel.Parameter(value=3.)
params['gamma_apf'] = radvel.Parameter(value=-3.7808)
params['jit_apf'] = radvel.Parameter(value=3.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0)
]
