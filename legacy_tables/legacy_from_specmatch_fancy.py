import pdb

import numpy as np
import pandas as pd
import radvel

specmatch = pd.read_csv('specmatch_results.csv') # Currently using Specmatch-synth
legacy = pd.read_csv('legacy_starlist.csv')

legacy_names = list(legacy['name'])

# Make specmatch names lowercase, to match the legacy target list file.
for i in np.arange(len(specmatch)):
    specmatch.loc[i, 'name'] = str(specmatch.loc[i, 'name']).lower()

#Make a new median-mass specmatch table.
specmatch_names = specmatch.name.unique()
specmatch_unique = pd.DataFrame(columns=['name', 'mass_median',
                                        'umstar1', 'umstar2', 'umstar'])

#Retrieve mass medians and uncertainties from specmatch table.
for i in np.arange(len(legacy_names)):
    legname = legacy_names[i]
    #spec = specmatch.query('name == {}'.format(legname))
    spec = specmatch.loc[specmatch['name'] == legname]
    mass_median = np.median(spec.iso_mass)
    umass1 = np.mean(spec.iso_mass_err1)
    umass2 = np.mean(spec.iso_mass_err2)
    mass_std = np.std(spec.iso_mass)
    specmatch_unique.loc[i] = [legname, mass_median, umass1, umass2, mass_std]

legacy_m = pd.merge(legacy, specmatch_unique, on=['name'])
legacy_m.to_csv('legacy_with_masses.csv')
