import pdb

import pandas as pd
import numpy as np
import radvel

import cpsutils
import rvsearch
from rvsearch import utils

starname = 'HD6872b'
nplanets = 1
instnames = ['k', 'j']
ntels = len(instnames)
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 2450000.

# load in data
#data = utils.read_from_csv('./setup_data/vst6872b.csv')
#data = data.query('jd < 2456000 or jd > 2457000')
data = cpsutils.io.loadcps('6872b', hires_rk=True, hires_rj=True, lick=False, ctslim=3000, binsize=0.5)
if 'jd' in data.columns:
    data['time'] = data['jd']
time_base = np.median(data['time'])
data['tel'] = data['tel'].str.decode('utf-8')

def initialize_params():
    params = radvel.Parameters(nplanets, basis='per tp e w k')

    params['per1'] = radvel.Parameter(value=15702.3)#11136.6)
    params['tp1'] = radvel.Parameter(value=2457218.85)#7135.72)
    params['k1'] = radvel.Parameter(value=5712.63)#7044.13)
    params['e1'] = radvel.Parameter(value=0.6605)#0.6936)
    params['w1'] = radvel.Parameter(value=-2.628)#-2.701)

    params['dvdt'] = radvel.Parameter(value=0, vary=False)
    params['curv'] = radvel.Parameter(value=0, vary=False)

    # Convert input orbital parameters into the fitting basis
    params = params.basis.to_any_basis(params,fitting_basis)

    return params

# initialize the orbit parameters and the orbit model
params = initialize_params()
params['gamma_j'] = radvel.Parameter(value=331.414, vary=False, linear=True)
params['jit_j'] = radvel.Parameter(value=1.)
params['gamma_k'] = radvel.Parameter(value=317.982, vary=False, linear=True)
params['jit_k'] = radvel.Parameter(value=2.)

priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_k', 0.0, 50.0),
    radvel.prior.HardBounds('jit_j', 0.0, 50.0)
]
