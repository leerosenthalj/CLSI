import radvel
import numpy as np
import pandas as pd

starname = 'HD 3167'
nplanets = 3
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c', 3: 'd'}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k',
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=0.959665)
anybasis_params['tc1'] = radvel.Parameter(value=2457744.651325)
anybasis_params['e1'] = radvel.Parameter(value=0)
anybasis_params['w1'] = radvel.Parameter(value=0.000000)
anybasis_params['k1'] = radvel.Parameter(value=3.580000)
anybasis_params['per2'] = radvel.Parameter(value=29.844982)
anybasis_params['tc2'] = radvel.Parameter(value=2457753.119284)
anybasis_params['e2'] = radvel.Parameter(value=0)
anybasis_params['w2'] = radvel.Parameter(value=0.000000)
anybasis_params['k2'] = radvel.Parameter(value=2.230000)
anybasis_params['per3'] = radvel.Parameter(value=8.509000)
anybasis_params['tc3'] = radvel.Parameter(value=2457743.990854)
anybasis_params['e3'] = radvel.Parameter(value=0.360000)
anybasis_params['w3'] = radvel.Parameter(value=0.000000)
anybasis_params['k3'] = radvel.Parameter(value=2.390000)

time_base = 2457743.990854
anybasis_params['dvdt'] = radvel.Parameter(value=0.0, vary=False)
anybasis_params['curv'] = radvel.Parameter(value=0.0, vary=False)
data = pd.read_csv('/data/radvel/input_dir/3167/data/all_data/all_data.csv',
                   dtype={'time': np.float64, 'mnvel': np.float64, 'err': np.float64, 'tel': str})
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.1)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel

instnames = ['apf', 'j']
ntels = len(instnames)
anybasis_params['gamma_apf'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_apf'] = radvel.Parameter(value=1.0)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=False, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)

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
mod.params['dvdt'].vary = False
mod.params['curv'].vary = False
mod.params['jit_apf'].vary = True
mod.params['jit_j'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          radvel.prior.Gaussian('per1', 0.959665, 2.1e-05),
          radvel.prior.Gaussian('per2', 29.844982, 0.001389),
          radvel.prior.Gaussian('per3', 8.509000, 0.045)
         ]
