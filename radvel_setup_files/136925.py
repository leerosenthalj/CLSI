import pdb

import pandas as pd
import numpy as np
import radvel

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

starname = 'HD136925'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
#data = utils.read_from_csv('setup_data/vst136925.csv')
data = cpsutils.io.loadcps('136925', hires_rk=True, hires_rj=True, lick=False, ctslim=3000, binsize=0.5)
data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    '''
    #Passing through the first datapoint.
    params = radvel.Parameters(1,basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=3900.74)
    params['tp1'] = radvel.Parameter(value=2458400.11) #8750 3048.11
    params['k1'] = radvel.Parameter(value=22.27)
    params['e1'] = radvel.Parameter(value=0.708) #0.7
    params['w1'] = radvel.Parameter(value=1.952) #-1.0
    params['dvdt'] = radvel.Parameter(value=0.0045, vary=True)
    params['curv'] = radvel.Parameter(value=1.107e-06, vary=True)
    '''
    #'''
    #BEST OPTION AT THE MOMENT
    #Not passing through the first datapoint. Try fitting proper periastron.
    params = radvel.Parameters(1,basis='per tc secosw sesinw k')
    params['per1'] = radvel.Parameter(value=4150)
    params['tc1'] = radvel.Parameter(value=2458750) #9000.11)#4450.11) #
    params['k1'] = radvel.Parameter(value=23.22)
    params['secosw1'] = radvel.Parameter(value=-0.33) #0.7
    params['sesinw1'] = radvel.Parameter(value=0.75) #-1.0
    params['dvdt'] = radvel.Parameter(value=0.0045, vary=True)
    params['curv'] = radvel.Parameter(value=1.107e-06, vary=True)
    params['dvdt'].vary = True
    params['curv'].vary = True
    #'''
    '''
    #Try best results from mcmc.
    params = radvel.Parameters(1,basis='per tc secosw sesinw k')
    params['per1'] = radvel.Parameter(value=4206.5)
    params['tc1'] = radvel.Parameter(value=2458400.)# 4300 is from mcmc, +per 8679.8
    params['k1'] = radvel.Parameter(value=22.4)
    params['secosw1'] = radvel.Parameter(value=-0.10) #0.7
    params['sesinw1'] = radvel.Parameter(value=0.83) #-1.0
    params['dvdt'] = radvel.Parameter(value=0.0045, vary=True)
    params['curv'] = radvel.Parameter(value=1.107e-06, vary=True)
    params['dvdt'].vary = True
    params['curv'].vary = True
    '''
    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)
    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-3.85, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=4.)
params['gamma_k'] = radvel.Parameter(value=2.71, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=4.)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]
