import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils.io

import rvsearch
from rvsearch import utils

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'HD156846'
nplanets = 1
instnames = ['j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.05, mstar_err=.06)

# load in data
'''
data =
'''
#data = utils.read_from_csv('setup_data/HIP52942A.csv')
data = cpsutils.io.loadcps('156846', hires_rk=False, hires_rj=True,
                           ctslim=3000, binsize=0.5)
data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=359.569)
    params['tc1'] = radvel.Parameter(value=2450763.63)
    params['k1'] = radvel.Parameter(value=462.977)
    params['e1'] = radvel.Parameter(value=0.48)
    params['w1'] = radvel.Parameter(value=0.89)

    # params['per2'] = radvel.Parameter(value=450.0)
    # params['tc2'] = radvel.Parameter(value=2450416.0)
    # params['k2'] = radvel.Parameter(value=3101.0)
    # params['e2'] = radvel.Parameter(value=0.8)
    # params['w2'] = radvel.Parameter(value=2.0)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)


    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-135.929, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=7.0, vary=False)


priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.PositiveKPrior( 1 ),
    radvel.prior.HardBounds('jit_j', 0.0, 30.0),
]
