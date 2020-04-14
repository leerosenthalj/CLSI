import radvel
import numpy as np
import pandas as pd

import cpsutils
import rvsearch
from rvsearch import utils

"""
"keywords"
"""
vary_dvdt = False # include a trend
"""
"""

starname = 'GJ 876'
nplanets = 2
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c', 3: 'd', 4: 'e'}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k',
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=61.116600)
anybasis_params['tc1'] = radvel.Parameter(value=2453056.9)
anybasis_params['e1'] = radvel.Parameter(value=0.043)
anybasis_params['w1'] = radvel.Parameter(value=-1.87)
anybasis_params['k1'] = radvel.Parameter(value=212.0)
anybasis_params['per2'] = radvel.Parameter(value=30.1934)
anybasis_params['tc2'] = radvel.Parameter(value=2453028.66)
anybasis_params['e2'] = radvel.Parameter(value=0.05)
anybasis_params['w2'] = radvel.Parameter(value=-1.22)
anybasis_params['k2'] = radvel.Parameter(value=88.340000)
#anybasis_params['per3'] = radvel.Parameter(value=1.937780)
#anybasis_params['tc3'] = radvel.Parameter(value=2453014.1)
#anybasis_params['e3'] = radvel.Parameter(value=0)
#anybasis_params['w3'] = radvel.Parameter(value=np.pi)
#anybasis_params['k3'] = radvel.Parameter(value=6.560000)
# anybasis_params['per4'] = radvel.Parameter(value=124.156)
# anybasis_params['tc4'] = radvel.Parameter(value=2453070.1)
# anybasis_params['e4'] = radvel.Parameter(value=0)
# anybasis_params['w4'] = radvel.Parameter(value=np.pi)
# anybasis_params['k4'] = radvel.Parameter(value=1.87)

time_base = 2453012.188702
anybasis_params['dvdt'] = radvel.Parameter(value=0.0)
anybasis_params['curv'] = radvel.Parameter(value=0.0)
'''
data = pd.read_csv('setup_data/gl876_all_data.csv',
                   dtype={'time': np.float64, 'mnvel': np.float64, 'err': np.float64, 'tel': str})
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.1)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel
'''
data = cpsutils.io.loadcps('gl876', hires_rk=True, hires_rj=True,
                           apf=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

instnames = ['j', 'k']
ntels = len(instnames)
anybasis_params['gamma_j'] = radvel.Parameter(value=-74.4, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=15.2)
anybasis_params['gamma_k'] = radvel.Parameter(value=-4.3, vary=False, linear=True)
anybasis_params['jit_k'] = radvel.Parameter(value=13.9)

params = anybasis_params.basis.to_any_basis(anybasis_params,fitting_basis)
mod = radvel.RVModel(params, time_base=time_base)

mod.params['per1'].vary = True
mod.params['tc1'].vary = True
mod.params['secosw1'].vary = True
mod.params['sesinw1'].vary = True
mod.params['per2'].vary = True
mod.params['tc2'].vary = True
mod.params['secosw2'].vary = True
mod.params['sesinw2'].vary = True
#mod.params['per3'].vary = True
#mod.params['tc3'].vary = True
#mod.params['secosw3'].vary = True
#mod.params['sesinw3'].vary = True
# mod.params['per4'].vary = True
# mod.params['tc4'].vary = True
# mod.params['secosw4'].vary = True
# mod.params['sesinw4'].vary = True
mod.params['dvdt'].vary = False
mod.params['curv'].vary = False
mod.params['jit_j'].vary = True
mod.params['jit_k'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          # radvel.prior.HardBounds('dvdt', -1.0, 1.0),
          # radvel.prior.Gaussian('per1', 61.116600, 0.0086),
          # radvel.prior.Gaussian('per2', 30.088100, 0.0082),
          # radvel.prior.Gaussian('per3', 1.937780, 2e-05),
          # radvel.prior.Gaussian('per4', 124.260000, 0.7),
          radvel.prior.HardBounds('jit_j', 0.0, 50.0),
          radvel.prior.HardBounds('jit_k', 0.0, 50.0)
         ]
