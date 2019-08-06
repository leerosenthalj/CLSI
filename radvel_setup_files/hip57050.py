import radvel
import numpy as np
import pandas as pd

starname = 'GJ 1148'
nplanets = 2
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'c', 2: 'b'}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k', 
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=532.580000)
anybasis_params['tc1'] = radvel.Parameter(value=2454847.762008)
anybasis_params['e1'] = radvel.Parameter(value=0.342000)
anybasis_params['w1'] = radvel.Parameter(value=3.672173)
anybasis_params['k1'] = radvel.Parameter(value=11.340000)
anybasis_params['per2'] = radvel.Parameter(value=41.397000)
anybasis_params['tc2'] = radvel.Parameter(value=2454847.762008)
anybasis_params['e2'] = radvel.Parameter(value=0.314000)
anybasis_params['w2'] = radvel.Parameter(value=4.155629)
anybasis_params['k2'] = radvel.Parameter(value=38.370000)

time_base = 2454847.762008
anybasis_params['dvdt'] = radvel.Parameter(value=0.0)
anybasis_params['curv'] = radvel.Parameter(value=0.0)
data = pd.read_csv('/data/radvel/input_dir/HIP57050/data/all_data/all_data.csv',
                   dtype={'time': np.float64, 'mnvel': np.float64, 'err': np.float64, 'tel': str})
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.1)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel

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
mod.params['per2'].vary = True
mod.params['tc2'].vary = True
mod.params['secosw2'].vary = True
mod.params['sesinw2'].vary = True
mod.params['dvdt'].vary = True
mod.params['curv'].vary = False
mod.params['jit_j'].vary = True
mod.params['jit_k'].vary = True

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          radvel.prior.Gaussian('per1', 532.580000, 2.52),
          radvel.prior.Gaussian('per2', 41.397000, 0.016)
         ]



