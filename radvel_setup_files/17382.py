import pdb

import pandas as pd
import numpy as np
import radvel
from cpsutils import io

import cpsutils
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

starname = 'HD17382'
nplanets = 1
instnames = ['j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=0.726, mstar_err=0.03)

# load in data
#data = utils.read_from_csv('./setup_data/vst17382.csv', binsize=0.5)
data = io.loadcps('17382', hires_rj=True, hires_rk=False,
                  verbose=False, ctslim=3000, detrend=False, binsize=1.0)
data['time'] = data['jd']
#data['tel'] = data['tel'].str.decode('utf-8')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(1, basis='per tc e w k')
    #From Halbwachs et al. 2018
    params['per1'] = radvel.Parameter(value=5554.71)
    params['tc1'] = radvel.Parameter(value=2453540.68)#2448024.2)
    params['k1'] = radvel.Parameter(value=2812.17)
    params['e1'] = radvel.Parameter(value=0.656)
    params['w1'] = radvel.Parameter(value=1.987)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params


# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-50.279, linear=True, vary=False)
params['jit_j'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
]
