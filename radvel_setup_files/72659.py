import radvel
import numpy as np
import pandas as pd
import cpsutils.io

import rvsearch
from rvsearch import utils

starname = 'HD 72659'
nplanets = 1
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b'}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k',
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=3383.000000)
anybasis_params['tc1'] = radvel.Parameter(value=2456919.019444)
anybasis_params['e1'] = radvel.Parameter(value=0)
anybasis_params['w1'] = radvel.Parameter(value=4.206243)
anybasis_params['k1'] = radvel.Parameter(value=42.400000)
anybasis_params['dvdt'] = radvel.Parameter(value=0, vary=True)
anybasis_params['curv'] = radvel.Parameter(value=0, vary=True)

data = cpsutils.io.loadcps('72659', hires_rk=True, hires_rj=True,
                           lick=True, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

instnames = ['j', 'k']
ntels = len(instnames)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_k'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_k'] = radvel.Parameter(value=1.0)

params = anybasis_params.basis.to_any_basis(anybasis_params,fitting_basis)
mod = radvel.RVModel(params, time_base=time_base)

mod.params['per1'].vary = True
mod.params['tc1'].vary = True
mod.params['secosw1'].vary = True
mod.params['sesinw1'].vary = True
mod.params['dvdt'].vary = True
mod.params['curv'].vary = True
mod.params['jit_j'].vary = True
mod.params['jit_k'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets)
         ]
