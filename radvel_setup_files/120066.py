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
linearP = True  # fit in linear P
informative_per_prior = False # include A Vandenburg's informative period prior
######
informative_per_prior_with_dur = False # same as above, but assume event duration is 2015-mid-2018 (3.5yrs)
fit_recentpoints_only = False # do a fit to only data taken after 2017
vary_dvdt = False # include a trend
"""
"""

starname = 'HD120066'
nplanets = 1
instnames = ['k', 'j', 'apf']#, 'mcdonald']
ntels = len(instnames)
fitting_basis = 'logper tc secosw sesinw logk'
if linearP:
    fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2440000.

# stellar mass & error
stellar = dict(mstar=1.16, mstar_err=.12)

# load in data
'''
data_folder = './merged_datasets/'
#data_folder = '/data/user/lrosenth/legacy/radvel_setup_files/merged_datasets/'
data = utils.read_from_csv(data_folder+'pi_with_mcd.csv', binsize=0.5)
data['time'] = data['jd']
'''
data_cps = cpsutils.io.loadcps('120066', hires_rk=True, hires_rj=True,
                               apf=True, ctslim=3000, binsize=0.5)
data['time'] = data['jd']
data['tel'] = data['tel'].str.decode('utf-8')
time_base = np.median(data.time.values)
'''
data_cps = cpsutils.io.loadcps('120066', hires_rk=True, hires_rj=True,
                               apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data_cps.columns:
    data_cps['time'] = data['jd']
data_mcd = pd.read_csv('/data/user/sblunt/Dropbox/120066_radvel/HD120066_McD.ALL',
	names=['time','mnvel','errvel', 'SVAL','sval_err'], header=None,sep='\s+'
)
data_mcd['tel']='m'
data_mcd['time'] -= 40000.
data = pd.concat([data_cps, data_mcd], ignore_index=True)
#Option to save merged dataset on cadence.
data.to_csv('/data/user/lrosenth/legacy/merged_datasets/120066_with_mcdonald.csv')
if fit_recentpoints_only:
    data = data[data.time >= 17754.]
    instnames = ['j', 'a', 'mcdonald']
'''
baseline = np.max(data.time.values) - np.min(data.time.values)


def initialize_params():
    # Sarah's period guess: 21860.
    params = radvel.Parameters(1,basis='per tc secosw sesinw k')
    # Guess #2
    params['per1'] = radvel.Parameter(value=21843.2)#21860.)
    params['tc1'] = radvel.Parameter(value=2458936.)#2458975.79)
    params['secosw1'] = radvel.Parameter(value=0.849)
    params['sesinw1'] = radvel.Parameter(value=-0.313)
    params['k1'] = radvel.Parameter(value=37.86)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)
    '''
    # Sarah's parameters
    params = radvel.Parameters(1,basis='per tp e w k')
    params['per1'] = radvel.Parameter(value=21860.)
    params['tp1'] = radvel.Parameter(value=2458104.64)#18134.)
    params['e1'] = radvel.Parameter(value=0.84)
    params['w1'] = radvel.Parameter(value=-0.26)
    params['k1'] = radvel.Parameter(value=38.)
    params['dvdt'] = radvel.Parameter(value=0, vary=vary_dvdt)
    params['curv'] = radvel.Parameter(value=0, vary=False)
    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)
    '''
    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=-44.57, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=2.56)
params['gamma_k'] = radvel.Parameter(value=-45.33, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=3.36)
params['gamma_apf'] = radvel.Parameter(value=-38.76, vary=False, linear=True)
params['jit_apf'] = radvel.Parameter(value=4.03)
#params['gamma_mcdonald'] = radvel.Parameter(value=-5.07, vary=False, linear=True)
#params['jit_mcdonald'] = radvel.Parameter(value=6.13)

priors = [
    radvel.prior.EccentricityPrior( 1 ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 10.0),
    radvel.prior.HardBounds('jit_j', 0.0, 10.0),
    radvel.prior.HardBounds('jit_apf', 0.0, 10.0)
]
#radvel.prior.HardBounds('jit_mcdonald', 0.0, 20.0)


if informative_per_prior:
   priors.append(radvel.prior.InformativeBaselinePrior('logper1', baseline))
elif informative_per_prior_with_dur:
    duration = 3.5*365.25
    priors.append(radvel.prior.InformativeBaselinePrior('logper1', baseline, duration))
