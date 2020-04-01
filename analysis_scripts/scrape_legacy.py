import os
import numpy as np
import pandas as pd
import radvel

from catalog_scrape import scrape

starnames = [name for name in os.listdir('.')
            if os.path.isdir(os.path.join('.', name))]

system_props = scrape(starnames, star_db_name='../CLSI/legacy_tables/legacy_specmatch_medians.csv')
