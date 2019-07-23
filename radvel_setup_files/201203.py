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

starname = 'HD 201203'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.4, mstar_err=.08)

# load in data
data = cpsutils.io.loadcps('201203', hires_rk=True, hires_rj=True, lick=False, ctslim=303, binsize=0.0)
#data = utils.read_from_csv('./setup_data/vst201203.csv', binsize=0.5)
data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=16469.7)
    params['tc1'] = radvel.Parameter(value=2457211.2)
    params['k1'] = radvel.Parameter(value=786.5)
    params['e1'] = radvel.Parameter(value=0.645)
    params['w1'] = radvel.Parameter(value=-2.48)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=5.0, vary=True)
params['gamma_k'] = radvel.Parameter(value=0.0, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=5.0, vary=True)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.1, 50.0),
    radvel.prior.HardBounds('jit_j', 0.1, 50.0),
]
