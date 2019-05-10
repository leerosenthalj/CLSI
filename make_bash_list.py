# Retrieve full list of legacy starnames, and write to
# bash list of search commands.
import os

import numpy as np
import pandas as pd

filename = 'legacy_analysis/legacy_starlist.csv'
starlist = pd.read_csv(filename)

bash_target = 'legacy_full_search.sh'
bash_writer = open(bash_target, 'w')

for name in starlist.name:
    bash_writer.write('python search_away.py \'{}\' \n'.format(name))
bash_writer.close()
