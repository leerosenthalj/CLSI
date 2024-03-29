import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD10790'
nplanets = 1
instnames = ['j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
#data = utils.read_from_csv('./setup_data/vst10790.csv')
data = cpsutils.io.loadcps('10790', hires_rk=False, hires_rj=True, lick=False, ctslim=3000, binsize=0.5)
#data = data[data['obnm'] != 'rj179.332']  # low counts, only 30k and also poor seeing
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=694.29)
    params['tp1'] = radvel.Parameter(value=2456800)
    params['k1'] = radvel.Parameter(value=8239.)
    params['e1'] = radvel.Parameter(value=0.84)
    params['w1'] = radvel.Parameter(value=-2.6)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=137, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]
