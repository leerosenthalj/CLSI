import os
import pickle
import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import radvel

# Draw end-of-chain samples from eccentricity posterior of each planet in each system in the catalog. Save as a table.
system_props = pd.read_csv('system_props_719.csv')

df_all = pd.DataFrame()
df_big = pd.DataFrame()
df_big_noncirc = pd.DataFrame()

df_big_multi = pd.DataFrame()
df_big_single = pd.DataFrame()

df_big_multi_noncirc = pd.DataFrame()
df_big_single_noncirc = pd.DataFrame()

df_cold = pd.DataFrame()
df_cold_small = pd.DataFrame()

df_giant_friend = pd.DataFrame()
df_giant_solo = pd.DataFrame()

stacker = []
names = []
dirs = [x[0] for x in os.walk('.')]
for dir in dirs:
    print(dir)
    name = dir[2:]
    pchain_path = dir + '/{}_pchains.csv'.format(name)
    if os.path.exists(pchain_path):
        system = system_props.query('name == "{}"'.format(name))
        print(name)
        names.append(name)
        chains = pd.read_csv(pchain_path)[-10000:]
        #chains = chains.sample(1000).reset_index()

        post = radvel.posterior.load(dir + '/post_final.pkl')
        nplanets = post.params.num_planets

        # Identify trends.
        keplerian = False
        for n in np.arange(1, nplanets+1):
            status = np.array(system['status{}'.format(n)])[0]
            if status == 'SS' or status == 'S':
                keplerian = True
        if post.params['dvdt'].value != 0. or keplerian:
            outer = True
        else:
            outer = False

        nplanets_true = 0
        nplanets_big = 0
        nplanets_big_noncirc = 0
        nplanets_hjs = 0

        for n in np.arange(1, nplanets+1):
            status = np.array(system['status{}'.format(n)])[0]
            mass = np.array(system['M{}'.format(n)])[0]
            axis = np.array(system['a{}'.format(n)])[0]

            if status == 'K' or status == 'C':
                nplanets_true += 1
                df_all['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                df_all['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                df_all['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                df_all['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                df_all['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

                if mass > 0.1:
                    nplanets_big += 1
                    df_big['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_big['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_big['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_big['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_big['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

                    if axis > 0.1:
                        nplanets_big_noncirc += 1
                        df_big_noncirc['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                        df_big_noncirc['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                        df_big_noncirc['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                        df_big_noncirc['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                        df_big_noncirc['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]
                    elif axis < 0.1:
                        nplanets_hjs += 1

                if mass > 0.3 and mass < 3 and axis > 1 and axis < 10:
                    df_cold['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_cold['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_cold['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_cold['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_cold['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]
                if mass > 0.1 and mass < 1 and axis > 1 and axis < 10:
                    df_cold_small['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_cold_small['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_cold_small['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_cold_small['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_cold_small['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]


        if nplanets_hjs == 0:
            for n in np.arange(1, nplanets+1):
                status = np.array(system['status{}'.format(n)])[0]
                mass = np.array(system['M{}'.format(n)])[0]
                axis = np.array(system['a{}'.format(n)])[0]
                good_status = (status == 'K' or status == 'C')

                if mass > 0.1 and axis > 0.1 and good_status:
                    df_giant_solo['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_giant_solo['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_giant_solo['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_giant_solo['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_giant_solo['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

        else:
            for n in np.arange(1, nplanets+1):
                status = np.array(system['status{}'.format(n)])[0]
                mass = np.array(system['M{}'.format(n)])[0]
                axis = np.array(system['a{}'.format(n)])[0]
                good_status = (status == 'K' or status == 'C')

                if mass > 0.1 and axis > 0.1 and good_status:
                    df_giant_friend['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_giant_friend['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_giant_friend['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_giant_friend['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_giant_friend['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]


        if nplanets_big > 1 or outer:
            for n in np.arange(1, nplanets+1):
                status = np.array(system['status{}'.format(n)])[0]
                mass = np.array(system['M{}'.format(n)])[0]
                axis = np.array(system['a{}'.format(n)])[0]
                good_status = (status == 'K' or status == 'C')

                if mass > 0.1 and good_status:
                    df_big_multi['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_big_multi['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_big_multi['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_big_multi['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_big_multi['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

                    if axis > 0.1:
                        df_big_multi_noncirc['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                        df_big_multi_noncirc['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                        df_big_multi_noncirc['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                        df_big_multi_noncirc['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                        df_big_multi_noncirc['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

        else:
            for n in np.arange(1, nplanets+1):
                status = np.array(system['status{}'.format(n)])[0]
                mass = np.array(system['M{}'.format(n)])[0]
                axis = np.array(system['a{}'.format(n)])[0]
                good_status = (status == 'K' or status == 'C')

                if mass > 0.1 and good_status:
                    df_big_single['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                    df_big_single['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                    df_big_single['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                    df_big_single['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                    df_big_single['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

                    if axis > 0.1:
                        df_big_single_noncirc['{0}_e{1}'.format(name, n)] = chains['e{}'.format(n)]
                        df_big_single_noncirc['{0}_a{1}'.format(name, n)] = chains['a{}'.format(n)]
                        df_big_single_noncirc['{0}_M{1}'.format(name, n)] = chains['M{}'.format(n)]
                        df_big_single_noncirc['{0}_k{1}'.format(name, n)] = chains['k{}'.format(n)]
                        df_big_single_noncirc['{0}_per{1}'.format(name, n)] = chains['per{}'.format(n)]

print('Done with bookkeeping')
df_all.to_csv('e_samples_all.csv')
print('Saved 1')
df_big.to_csv('e_samples_big.csv')
df_big_noncirc.to_csv('e_samples_big_noncirc.csv')
print('Saved 3')
df_big_multi.to_csv('e_samples_big_multi.csv')
df_big_single.to_csv('e_samples_big_single.csv')

df_big_multi_noncirc.to_csv('e_samples_big_multi_noncirc.csv')
df_big_single_noncirc.to_csv('e_samples_big_single_noncirc.csv')

df_cold.to_csv('e_samples_cold.csv')
df_cold_small.to_csv('e_samples_cold_small.csv')

df_giant_friend.to_csv('e_samples_giant_friend.csv')
df_giant_solo.to_csv('e_samples_giant_solo.csv')
