import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD1461'
nplanets = 2
instnames = ['j', 'apf']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
data = utils.read_from_csv('../../rvdata/vst1461.csv')
data = data.query('tel != "k"').reset_index()
#data = cpsutils.io.loadcps('156668', hires_rk=False, hires_rj=True,
#                           apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
#data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=5.77124)
    params['tp1'] = radvel.Parameter(value=2455722.05)
    params['k1'] = radvel.Parameter(value=2.25922)
    params['e1'] = radvel.Parameter(value=0.205714)
    params['w1'] = radvel.Parameter(value=-0.441632)
    '''
    params['per2'] = radvel.Parameter(value=13.5053)
    params['tp2'] = radvel.Parameter(value=2455204.82)
    params['k2'] = radvel.Parameter(value=2.19857)
    params['e2'] = radvel.Parameter(value=0.146231)
    params['w2'] = radvel.Parameter(value=-2.43268)
    '''
    # ACTIVITY, NEED TO TEST
    params['per2'] = radvel.Parameter(value=3482.66)
    params['tp2'] = radvel.Parameter(value=2456818.23)
    params['k2'] = radvel.Parameter(value=3.0)
    params['e2'] = radvel.Parameter(value=0.3983)
    params['w2'] = radvel.Parameter(value=-0.494308)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=3.904, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=3.)
#params['gamma_k'] = radvel.Parameter(value=3.904, vary=False, linear=True)
#params['jit_k'] = radvel.Parameter(value=3.)
params['gamma_apf'] = radvel.Parameter(value=-3.7808, vary=False, linear=True)
params['jit_apf'] = radvel.Parameter(value=3.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0)
]
#    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
