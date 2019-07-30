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

starname = 'HD16287'
nplanets = 1
instnames = ['j', 'k']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
# stellar = dict(mstar=1.05, mstar_err=.06)

# load in data
'''
data =
'''
data = cpsutils.io.loadcps('16287', hires_rk=True, hires_rj=True,
                           ctslim=4000, binsize=0.5)
data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=14.8386, vary=True)
    params['tc1'] = radvel.Parameter(value=2450770.578)
    params['k1'] = radvel.Parameter(value=10713.9)
    params['e1'] = radvel.Parameter(value=0.206929)
    params['w1'] = radvel.Parameter(value=0.179116)

    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)


    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-1541.31, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=7.0, vary=False)
params['gamma_k'] = radvel.Parameter(value=-1534.38, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=7.0, vary=False)


priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.PositiveKPrior( 1 ),
]
