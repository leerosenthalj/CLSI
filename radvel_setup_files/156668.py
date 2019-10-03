import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD156668'
nplanets = 3
instnames = ['j', 'apf']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('../../rvdata/vst156668.csv')
data = data.query('tel != "k"')
#data = cpsutils.io.loadcps('156668', hires_rk=False, hires_rj=True,
#                           apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
#data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=4.6428)
    params['tp1'] = radvel.Parameter(value=2455202.95)
    params['k1'] = radvel.Parameter(value=2.475)
    params['e1'] = radvel.Parameter(value=0.281845)
    params['w1'] = radvel.Parameter(value=-2.490)

    params['per2'] = radvel.Parameter(value=811.459)
    params['tp2'] = radvel.Parameter(value=2457743.68)
    params['k2'] = radvel.Parameter(value=2.61)
    params['e2'] = radvel.Parameter(value=0.161)
    params['w2'] = radvel.Parameter(value=-2.846)

    # ACTIVITY, NEED TO TEST
    #params['per3'] = radvel.Parameter(value=3470.0)
    params['per3'] = radvel.Parameter(value=3470.0)
    params['tp3'] = radvel.Parameter(value=2458556.72)
    params['k3'] = radvel.Parameter(value=3.05)
    params['e3'] = radvel.Parameter(value=0.25)
    params['w3'] = radvel.Parameter(value=-2.23)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=3.904, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=3.)
params['gamma_apf'] = radvel.Parameter(value=-3.7808, vary=False, linear=True)
params['jit_apf'] = radvel.Parameter(value=3.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0)
]
