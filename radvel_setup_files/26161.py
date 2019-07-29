import pandas as pd
import numpy as np
import radvel
import cpsutils

from rvsearch import utils

# radvel fit -s 120066.py -d 2019-2-19/master
# radvel mcmc -s 120066.py -d 2019-2-19/master --maxGR 1.001 --minsteps 1000 --nsteps 10000 --minpercent 100

"""
"keywords"
"""
linearP = True  # fit in linear P
fit_recentpoints_only = False

starname = 'HD26161'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
#fitting_basis = 'logper tc secosw sesinw logk'
#if linearP:
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2440000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
#data = utils.read_from_csv('setup_data/vst26161.csv')
data = cpsutils.io.loadcps('26161', hires_rk=True, hires_rj=True,
                           apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

baseline = np.max(data.time.values) - np.min(data.time.values)

def initialize_params():
    params = radvel.Parameters(1,basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=18930.7)
    params['tp1'] = radvel.Parameter(value=2458730.1)
    params['e1'] = radvel.Parameter(value=0.74)
    params['w1'] = radvel.Parameter(value=0.326)
    params['k1'] = radvel.Parameter(value=47.63)
    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-40., vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=2.)
params['gamma_k'] = radvel.Parameter(value=-40., vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=3.)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]
