import pdb

import numpy as np
import pandas as pd
import radvel

import rvsearch
from rvsearch import utils

specmatch = pd.read_csv('specmatch_results.csv')
legacy = pd.read_csv('legacy_starlist.csv')

legacy_names = list(legacy['name'])

# Make specmatch names lowercase, to match the legacy target list file.
for i in np.arange(len(specmatch)):
    specmatch.loc[i, 'name'] = str(specmatch.loc[i, 'name']).lower()

#Make a new median-mass specmatch table.
specmatch_names = specmatch.name.unique()
specmatch_unique = pd.DataFrame(columns=['name', 'mass_mean', 'mass_median',
                                        'umass1', 'umass2', 'mass_std'])
for i in np.arange(len(specmatch_names)):
    name = specmatch_names[i]
    spec_masses = specmatch.loc[specmatch['name'] == name]['iso_mass']
    spec_umass1 = specmatch.loc[specmatch['name'] == name]['iso_mass_err1']
    spec_umass2 = specmatch.loc[specmatch['name'] == name]['iso_mass_err2']
    mass_mean = np.mean(spec_masses)
    mass_median = np.median(spec_masses)
    umass1 = np.mean(spec_umass1)
    umass2 = np.mean(spec_umass2)
    mass_std = np.std(spec_masses)

    specmatch_unique.loc[i] = [name, mass_mean, mass_median, umass1, umass2, mass_std]

specmatch_unique.to_csv('specmatch_unique_props.csv')

legacy['mstar'] = np.nan
legacy['umstar'] = np.nan
legacy['umstar1'] = np.nan
legacy['umstar2'] = np.nan

#specmatch_unique = pd.read_csv('specmatch_unique_masses.csv')
for name in legacy_names:
    #pdb.set_trace()
    leg = legacy.index[legacy['name'] == str(name)][0]
    spec = specmatch_unique.index[specmatch_unique['name'] == str(name)]#[0]
    legacy.loc[leg, 'mstar'] = specmatch_unique.loc[spec, 'mass_median']
    legacy.loc[leg, 'umstar'] = specmatch_unique.loc[spec, 'mass_std']
    legacy.loc[leg, 'umstar1'] = specmatch_unique.loc[spec, 'umass1']
    legacy.loc[leg, 'umstar2'] = specmatch_unique.loc[spec, 'umass2']

legacy.to_csv('legacy_specmatch.csv')
