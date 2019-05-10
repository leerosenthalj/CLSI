import pdb

import pandas as pd
import numpy as np
import radvel

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

starname = 'HD4747'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
'''
data_cps = pd.read_csv('~/Dropbox/120066_radvel/120066.txt')
data_mcd = pd.read_csv('~/Dropbox/120066_radvel/HD120066_McD.ALL',
	names=['time','mnvel','errvel', 'SVAL','sval_err'], header=None,sep='\s+'
)
data_mcd['tel']='m'
data_mcd['time'] -= 40000.
'''
#data = pd.concat([data_cps, data_mcd], ignore_index=True)
data = utils.read_from_csv('vst4747.csv')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(1,basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=12525.0)
    params['tp1'] = radvel.Parameter(value=2450593.5)
    params['k1'] = radvel.Parameter(value=725.3)
    params['e1'] = radvel.Parameter(value=0.74)
    params['w1'] = radvel.Parameter(value=4.6967)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params
    #params['per1'] = radvel.Parameter(value=13835.67)
    #params['tp1'] = radvel.Parameter(value=2450593.5)
	#params['k1'] = radvel.Parameter(value=755.3)

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-214.9)
params['jit_j'] = radvel.Parameter(value=4.)
params['gamma_k'] = radvel.Parameter(value=-201.4)
params['jit_k'] = radvel.Parameter(value=4.)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
]
