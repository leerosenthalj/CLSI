import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD68988'
nplanets = 2
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
#data = utils.read_from_csv('./setup_data/vst68988.csv')
data = cpsutils.io.loadcps('68988', hires_rk=True, hires_rj=True,
                           ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=6.27642)
    params['tp1'] = radvel.Parameter(value=2451548.84)
    params['k1'] = radvel.Parameter(value=193.343)
    params['e1'] = radvel.Parameter(value=0.152965)
    params['w1'] = radvel.Parameter(value=0.67475)

    params['per2'] = radvel.Parameter(value=15408.5)
    params['tp2'] = radvel.Parameter(value=2453380.75)
    params['k2'] = radvel.Parameter(value=125.813)
    params['e2'] = radvel.Parameter(value=0.401936)
    params['w2'] = radvel.Parameter(value=3.08)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=86.8118, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=0.2)
params['gamma_k'] = radvel.Parameter(value=97.5783, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=1.0)
#params['gamma_lick'] = radvel.Parameter(value=222.728)
#params['jit_lick'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 2.0)
    #radvel.prior.HardBounds('jit_lick', 0.0, 10.0)
]
