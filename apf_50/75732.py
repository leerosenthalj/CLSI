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

starname = '55Cnc'
nplanets = 6
instnames = ['k', 'j', 'apf', 'lick']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data, define time-base.
#data = pd.concat([data_cps, data_mcd], ignore_index=True)
data = utils.read_from_csv('vst75732.csv')
time_base = np.median(data['time'])

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tc e w k')

    params['per1'] = radvel.Parameter(value=0.73654737, vary=False)
    params['tc1'] = radvel.Parameter(value=55733.006)
    params['e1'] = radvel.Parameter(value=0.05)
    params['w1'] = radvel.Parameter(value=1.501)
    params['k1'] = radvel.Parameter(value=6.02)

    params['per2'] = radvel.Parameter(value=14.6516, vary=False)
    params['tc2'] = radvel.Parameter(value=55495.587)
    params['e2'] = radvel.Parameter(value=0.001)
    params['w2'] = radvel.Parameter(value=-0.37525)
    params['k2'] = radvel.Parameter(value=71.37)

    params['per3'] = radvel.Parameter(value=44.3989, vary=False)
    params['tc3'] = radvel.Parameter(value=55492.02)
    params['e3'] = radvel.Parameter(value=0.03)
    params['w3'] = radvel.Parameter(value=0.04189)
    params['k3'] = radvel.Parameter(value=9.89)

    params['per4'] = radvel.Parameter(value=259.88, vary=False)
    params['tc4'] = radvel.Parameter(value=55491.5)
    params['e4'] = radvel.Parameter(value=0.08)
    params['w4'] = radvel.Parameter(value=-1.703)
    params['k4'] = radvel.Parameter(value=5.14)

    # Activity
    params['per5'] = radvel.Parameter(value=3822.4, vary=False)
    params['tc5'] = radvel.Parameter(value=553366.9)
    params['e5'] = radvel.Parameter(value=0.17)
    params['w5'] = radvel.Parameter(value=3.049)
    params['k5'] = radvel.Parameter(value=15.2)

    params['per6'] = radvel.Parameter(value=5574.2, vary=False)
    params['tc6'] = radvel.Parameter(value=56669.3)
    params['e6'] = radvel.Parameter(value=0.13)
    params['w6'] = radvel.Parameter(value=-1.206)
    params['k6'] = radvel.Parameter(value=38.6)

    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-23.6937)
params['jit_j'] = radvel.Parameter(value=5.)
params['gamma_k'] = radvel.Parameter(value=-29.2488)
params['jit_k'] = radvel.Parameter(value=5.)
params['gamma_apf'] = radvel.Parameter(value=9.09471)
params['jit_apf'] = radvel.Parameter(value=5.)
params['gamma_lick'] = radvel.Parameter(value=-0.170913)
params['jit_lick'] = radvel.Parameter(value=5.)

priors = [
    radvel.prior.EccentricityPrior(nplanets), # Keeps eccentricity < 1
	radvel.prior.PositiveKPrior(nplanets),
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0),
    radvel.prior.HardBounds('jit_lick', 0.0, 20.0)
]
