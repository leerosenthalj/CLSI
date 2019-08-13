import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD159868'
nplanets = 2
instnames = ['j', 'AAT']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('./merged_datasets/159868_keck_aat.csv')
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
#data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=1191.71)
    params['tp1'] = radvel.Parameter(value=2453913.38)
    params['k1'] = radvel.Parameter(value=38.65)
    params['e1'] = radvel.Parameter(value=0.0543)
    params['w1'] = radvel.Parameter(value=2.42)

    params['per2'] = radvel.Parameter(value=350.164)
    params['tp2'] = radvel.Parameter(value=2454989.04)
    params['k2'] = radvel.Parameter(value=20.201)
    params['e2'] = radvel.Parameter(value=0.154)
    params['w2'] = radvel.Parameter(value=-1.38)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-1.90, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=4.2)
params['gamma_AAT'] = radvel.Parameter(value=6.87, vary=False, linear=True)
params['jit_AAT'] = radvel.Parameter(value=7.37)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_AAT', 0.0, 10.0)
]
