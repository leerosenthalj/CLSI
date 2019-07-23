import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD26990'
nplanets = 1
instnames = ['j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
#data = utils.read_from_csv('./setup_data/vst26990.csv')
data = cpsutils.io.loadcps('26990', hires_rk=False, hires_rj=True,
                           lick=False, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=2076.46)
    params['tp1'] = radvel.Parameter(value=2454625.2)
    params['k1'] = radvel.Parameter(value=8018.54)
    params['e1'] = radvel.Parameter(value=0.83)
    params['w1'] = radvel.Parameter(value=0.73)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-2321.)
params['jit_j'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]
