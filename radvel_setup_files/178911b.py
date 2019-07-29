import pdb

import pandas as pd
import numpy as np
import radvel
from cpsutils import io

import rvsearch
from rvsearch import utils

# radvel fit -s 120066.py -d 2019-2-19/master
# radvel mcmc -s 120066.py -d 2019-2-19/master --maxGR 1.001 --minsteps 1000 --nsteps 10000 --minpercent 100

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'HD178911b'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=0.726, mstar_err=0.03)

# load in data
data = io.loadcps('178911b', hires_rj=True, hires_rk=True,
                  verbose=False, ctslim=3000, detrend=False, binsize=0.5)
data['time'] = data['jd']
data['tel'] = data['tel'].str.decode('utf-8')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(1, basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=71.4808)
    params['tp1'] = radvel.Parameter(value=2450300.608)
    params['k1'] = radvel.Parameter(value=343.232)
    params['e1'] = radvel.Parameter(value=0.11316)
    params['w1'] = radvel.Parameter(value=2.89598)
    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params


# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=57.8914, linear=True, vary=False)
params['jit_j'] = radvel.Parameter(value=2.)
params['gamma_k'] = radvel.Parameter(value=57.6918, linear=True, vary=False)
params['jit_k'] = radvel.Parameter(value=2.)


priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_k', 0.0, 10.0)
]
