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
vary_dvdt = True # include a trend
"""
"""

starname = 'GL 570B'
nplanets = 1
instnames = ['lick', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=0.57, mstar_err=.02)

# load in data
'''
data_cps = pd.read_csv('~/Dropbox/120066_radvel/120066.txt')
data_mcd = pd.read_csv('~/Dropbox/120066_radvel/HD120066_McD.ALL',
	names=['time','mnvel','errvel', 'SVAL','sval_err'], header=None,sep='\s+'
)
data_mcd['tel']='m'
data_mcd['time'] -= 40000.
'''
data = cpsutils.io.loadcps('gl570b', hires_rk=False, hires_rj=True, lick=True, ctslim=500, binsize=0.0)
data = data[data['obnm'] != 'rj20.427']  # these appear to be outliers
data = data[data['obnm'] != 'rj21.410']  # these appear to be outliers

data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(1,basis='per tc e w k')
    params['per1'] = radvel.Parameter(value=323.33)
    params['tc1'] = radvel.Parameter(value=2450305.3)
    params['k1'] = radvel.Parameter(value=27019.0)
    params['e1'] = radvel.Parameter(value=0.832)
    params['w1'] = radvel.Parameter(value=2.2)
    params['dvdt'] = radvel.Parameter(value=4.07, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-3963.8, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=1639.0, vary=True)
params['gamma_lick'] = radvel.Parameter(value=15698.8, vary=False, linear=True)
params['jit_lick'] = radvel.Parameter(value=959.0, vary=True)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_lick', 50.0, 15000.0),
#    radvel.prior.HardBounds('jit_j', 0.0, 30.0),
]
