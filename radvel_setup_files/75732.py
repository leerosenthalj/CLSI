import radvel
import numpy as np
import pandas as pd
from rvsearch import utils
import cpsutils.io

starname = '75732'
nplanets = 5
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', }

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k',
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=14.648000)
anybasis_params['tc1'] = radvel.Parameter(value=2455805.465756)
anybasis_params['e1'] = radvel.Parameter(value=0)
anybasis_params['w1'] = radvel.Parameter(value=0.715585)
anybasis_params['k1'] = radvel.Parameter(value=77.100000)
anybasis_params['per2'] = radvel.Parameter(value=44.380000)
anybasis_params['tc2'] = radvel.Parameter(value=2455846.148111)
anybasis_params['e2'] = radvel.Parameter(value=0)
anybasis_params['w2'] = radvel.Parameter(value=6.213372)
anybasis_params['k2'] = radvel.Parameter(value=10.120000)
anybasis_params['per3'] = radvel.Parameter(value=4909.000000)
anybasis_params['tc3'] = radvel.Parameter(value=2456162.677778)
anybasis_params['e3'] = radvel.Parameter(value=0)
anybasis_params['w3'] = radvel.Parameter(value=4.433136)
anybasis_params['k3'] = radvel.Parameter(value=45.200000)
anybasis_params['per4'] = radvel.Parameter(value=0.737000)
anybasis_params['tc4'] = radvel.Parameter(value=2455802.360600)
anybasis_params['e4'] = radvel.Parameter(value=0)
anybasis_params['w4'] = radvel.Parameter(value=1.570796)
anybasis_params['k4'] = radvel.Parameter(value=6.000000)
anybasis_params['per5'] = radvel.Parameter(value=261.200000)
anybasis_params['tc5'] = radvel.Parameter(value=2456021.647778)
anybasis_params['e5'] = radvel.Parameter(value=0.320000)
anybasis_params['w5'] = radvel.Parameter(value=2.426008)
anybasis_params['k5'] = radvel.Parameter(value=6.200000)

anybasis_params['dvdt'] = radvel.Parameter(value=0.0)
anybasis_params['curv'] = radvel.Parameter(value=0.0)

data = cpsutils.io.loadcps('75732', hires_rk=True, hires_rj=True, lick=True,
                           ctslim=3000, binsize=0.5)
data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

instnames = ['apf', 'j', 'k', 'lick']
ntels = len(instnames)
anybasis_params['gamma_apf'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_apf'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_k'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_k'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_lick'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_lick'] = radvel.Parameter(value=1.0)

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
mod.params['per3'].vary = True
mod.params['tc3'].vary = True
mod.params['secosw3'].vary = True
mod.params['sesinw3'].vary = True
mod.params['per4'].vary = True
mod.params['tc4'].vary = True
mod.params['secosw4'].vary = False
mod.params['sesinw4'].vary = False
mod.params['per5'].vary = True
mod.params['tc5'].vary = True
mod.params['secosw5'].vary = False
mod.params['sesinw5'].vary = False
mod.params['dvdt'].vary = False
mod.params['curv'].vary = False
mod.params['jit_apf'].vary = True
mod.params['jit_j'].vary = True
mod.params['jit_k'].vary = True
mod.params['jit_lick'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets)]#,
          #radvel.prior.Gaussian('per1', 14.648000, 0.0009),
          #radvel.prior.Gaussian('per2', 44.380000, 0.007),
          #radvel.prior.Gaussian('per3', 4909.000000, 30),
          #radvel.prior.Gaussian('per4', 0.737000, 3e-06),
          #radvel.prior.Gaussian('per5', 261.200000, 0.4),
          #radvel.prior.UserDefinedPrior(['gamma_j', 'gamma_k'], utils.GaussianDiffFunc, 'Gaussian Prior on HIRES offset')]

stellar = dict(mstar=0.9859, mstar_err=0.0405)
