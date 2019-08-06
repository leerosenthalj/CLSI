import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD75732'
nplanets = 5
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f'}

# Define prior centers (initial guesses) in a basis of your choice.
# initialize Parameters object
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k',
                                    planet_letters=planet_letters)
anybasis_params['per1'] = radvel.Parameter(value=14.651520)
anybasis_params['tc1'] = radvel.Parameter(value=2456340.313479)
anybasis_params['e1'] = radvel.Parameter(value=0.003400)
anybasis_params['w1'] = radvel.Parameter(value=1.710423)
anybasis_params['k1'] = radvel.Parameter(value=71.400000)
anybasis_params['per2'] = radvel.Parameter(value=44.417500)
anybasis_params['tc2'] = radvel.Parameter(value=2456340.313479)
anybasis_params['e2'] = radvel.Parameter(value=0.020000)
anybasis_params['w2'] = radvel.Parameter(value=0.890118)
anybasis_params['k2'] = radvel.Parameter(value=10.180000)
anybasis_params['per3'] = radvel.Parameter(value=4825.000000)
anybasis_params['tc3'] = radvel.Parameter(value=2456340.313479)
anybasis_params['e3'] = radvel.Parameter(value=0.019000)
anybasis_params['w3'] = radvel.Parameter(value=0.767945)
anybasis_params['k3'] = radvel.Parameter(value=48.290000)
anybasis_params['per4'] = radvel.Parameter(value=0.736539)
anybasis_params['tc4'] = radvel.Parameter(value=2456340.657675)
anybasis_params['e4'] = radvel.Parameter(value=0.000000)
anybasis_params['w4'] = radvel.Parameter(value=0.000000)
anybasis_params['k4'] = radvel.Parameter(value=0.000000)
anybasis_params['per5'] = radvel.Parameter(value=262.000000)
anybasis_params['tc5'] = radvel.Parameter(value=2456340.313479)
anybasis_params['e5'] = radvel.Parameter(value=0.305000)
anybasis_params['w5'] = radvel.Parameter(value=2.897247)
anybasis_params['k5'] = radvel.Parameter(value=4.870000)

time_base = 2456340.313479
anybasis_params['dvdt'] = radvel.Parameter(value=0.0)
anybasis_params['curv'] = radvel.Parameter(value=0.0)
data = pd.read_csv('/data/radvel/input_dir/75732/data/dynamic_custom/remove_outliers.csv')
print(data)
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.1)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel

instnames = ['apf', 'j', 'lick']
ntels = len(instnames)
anybasis_params['gamma_apf'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_apf'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_lick'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_lick'] = radvel.Parameter(value=1.0)

anybasis_params['dvdt'] = radvel.Parameter(value=0, vary=False)
anybasis_params['curv'] = radvel.Parameter(value=0, vary=False)

params = anybasis_params.basis.to_any_basis(anybasis_params,fitting_basis)
mod = radvel.RVModel(params, time_base=time_base)
'''
mod.params['per1'].vary = True
mod.params['tc1'].vary = True
mod.params['secosw1'].vary = True
mod.params['sesinw1'].vary = True
mod.params['per2'].vary = True
mod.params['tc2'].vary = True
mod.params['secosw2'].vary = True
mod.params['sesinw2'].vary = True
mod.params['per3'].vary = True
mod.params['tc3'].vary = True
mod.params['secosw3'].vary = True
mod.params['sesinw3'].vary = True
mod.params['per4'].vary = False
mod.params['tc4'].vary = False
mod.params['secosw4'].vary = True
mod.params['sesinw4'].vary = True
mod.params['per5'].vary = True
mod.params['tc5'].vary = True
mod.params['secosw5'].vary = True
mod.params['sesinw5'].vary = True
mod.params['dvdt'].vary = False
mod.params['curv'].vary = False
mod.params['gamma_apf'].vary = True
mod.params['jit_apf'].vary = True
mod.params['gamma_j'].vary = True
mod.params['jit_j'].vary = True
'''
priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets)
         ]
