import radvel
import numpy as np
import pandas as pd
import cpsutils.io

starname = 'HD 187123'
nplanets = 2
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c'}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k', 
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=3.096589)
anybasis_params['tc1'] = radvel.Parameter(value=2455665.887479)
anybasis_params['e1'] = radvel.Parameter(value=0)
anybasis_params['w1'] = radvel.Parameter(value=6.283185)
anybasis_params['k1'] = radvel.Parameter(value=68.910000)
anybasis_params['per2'] = radvel.Parameter(value=3324.000000)
anybasis_params['tc2'] = radvel.Parameter(value=2458717.183333)
anybasis_params['e2'] = radvel.Parameter(value=0)
anybasis_params['w2'] = radvel.Parameter(value=4.511676)
anybasis_params['k2'] = radvel.Parameter(value=25.100000)

time_base = 2455663.735766
anybasis_params['dvdt'] = radvel.Parameter(value=-0.00124)
anybasis_params['curv'] = radvel.Parameter(value=0.0)
data = cpsutils.io.loadcps('187123', hires_rk=True, hires_rj=True, apf=True, ctslim=3000, binsize=2.0)
data['time'] = data['jd']
data['tel'] = data['tel'].str.decode('utf-8')

instnames = ['apf', 'j', 'k']
ntels = len(instnames)
anybasis_params['gamma_apf'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_apf'] = radvel.Parameter(value=2.1)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=2.7)
anybasis_params['gamma_k'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_k'] = radvel.Parameter(value=2.2)

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
mod.params['dvdt'].vary = True
mod.params['curv'].vary = False
mod.params['jit_apf'].vary = True
mod.params['jit_j'].vary = True
mod.params['jit_k'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          # radvel.prior.HardBounds('jit_apf', 0.0, 50.0),
          # radvel.prior.HardBounds('jit_j', 0.0, 50.0),
          # radvel.prior.HardBounds('jit_k', 0.0, 50.0)
         ]



