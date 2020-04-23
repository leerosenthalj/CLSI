# Make a scatter plot of planet mass and semi-major axis.
import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import radvel

def make_trend_list(props_filename, save=True):
    props = pd.read_csv(props_filename)

    num_stars = len(props)
    starnames  = props['name']
    mstar = np.asarray(props.Mstar)
    dvdt       = np.asarray(props.dvdt)
    dvdt_med   = np.asarray(props.dvdt_med)
    dvdt_minus = np.asarray(props.dvdt_minus)
    dvdt_plus  = np.asarray(props.dvdt_plus)
    curv       = np.asarray(props.curv)
    curv_med   = np.asarray(props.curv_med)
    curv_minus = np.asarray(props.curv_minus)
    curv_plus  = np.asarray(props.curv_plus)

    trends_dict = {'hostname':starnames, 'mstar':mstar,'dvdt':dvdt,
                   'dvdt_med':dvdt_med, 'dvdt_minus':dvdt_minus,
                   'dvdt_plus':dvdt_plus, 'curv':curv, 'curv_med':curv_med,
                   'curv_minus':curv_minus,'curv_plus':curv_plus}
    trends = pd.DataFrame(trends_dict)
    trends = trends.loc[~np.isnan(trends.dvdt)].reset_index()
    trends.to_csv('trend_list.csv')

def make_planet_list(props_filename, save=True):
    props = pd.read_csv(props_filename)

    # Collect the starname/mass/axis sets for each planet.
    num_stars = len(props)
    starnames    = []
    masses       = []
    masses_med   = []
    masses_minus = []
    masses_plus  = []
    ks           = []
    ks_med       = []
    ks_minus     = []
    ks_plus      = []
    axes         = []
    axes_med     = []
    axes_minus   = []
    axes_plus    = []
    pers         = []
    pers_med     = []
    pers_minus   = []
    pers_plus    = []
    tcs          = []
    tcs_med      = []
    tcs_minus    = []
    tcs_plus     = []
    tps          = []
    tps_med      = []
    tps_minus    = []
    tps_plus     = []
    es           = []
    es_med       = []
    es_minus     = []
    es_plus      = []
    es_mode      = []
    es_68        = []
    ws           = []
    ws_med       = []
    ws_minus     = []
    ws_plus      = []
    insols       = []
    insols_med   = []
    insols_minus = []
    insols_plus  = []
    teqs         = []
    teqs_med     = []
    teqs_minus   = []
    teqs_plus    = []
    post_paths   = []

    for i in np.arange(num_stars):
        num_planets = props.loc[i, 'num_planets']
        if num_planets != 0:
            for j in np.arange(1, num_planets+1):
                starname = props.loc[i, 'name']
                starnames.append(starname)

                masses.append(props.loc[i, 'M{}'.format(j)])
                masses_med.append(props.loc[i, 'M{}_med'.format(j)])
                masses_minus.append(props.loc[i, 'M{}_minus'.format(j)])
                masses_plus.append(props.loc[i, 'M{}_plus'.format(j)])

                ks.append(props.loc[i, 'k{}'.format(j)])
                ks_med.append(props.loc[i, 'k{}_med'.format(j)])
                ks_minus.append(props.loc[i, 'k{}_minus'.format(j)])
                ks_plus.append(props.loc[i, 'k{}_plus'.format(j)])

                axes.append(props.loc[i, 'a{}'.format(j)])
                axes_med.append(props.loc[i, 'a{}_med'.format(j)])
                axes_minus.append(props.loc[i, 'a{}_minus'.format(j)])
                axes_plus.append(props.loc[i, 'a{}_plus'.format(j)])

                pers.append(props.loc[i, 'per{}'.format(j)])
                pers_med.append(props.loc[i, 'per{}_med'.format(j)])
                pers_minus.append(props.loc[i, 'per{}_minus'.format(j)])
                pers_plus.append(props.loc[i, 'per{}_plus'.format(j)])

                tcs.append(props.loc[i, 'tc{}'.format(j)])
                tcs_med.append(props.loc[i, 'tc{}_med'.format(j)])
                tcs_minus.append(props.loc[i, 'tc{}_minus'.format(j)])
                tcs_plus.append(props.loc[i, 'tc{}_plus'.format(j)])

                tps.append(props.loc[i, 'tp{}'.format(j)])
                tps_med.append(props.loc[i, 'tp{}_med'.format(j)])
                tps_minus.append(props.loc[i, 'tp{}_minus'.format(j)])
                tps_plus.append(props.loc[i, 'tp{}_plus'.format(j)])

                ws.append(props.loc[i, 'w{}'.format(j)])
                ws_med.append(props.loc[i, 'w{}_med'.format(j)])
                ws_minus.append(props.loc[i, 'w{}_minus'.format(j)])
                ws_plus.append(props.loc[i, 'w{}_plus'.format(j)])

                es.append(props.loc[i, 'e{}'.format(j)])
                es_med.append(props.loc[i, 'e{}_med'.format(j)])
                es_minus.append(props.loc[i, 'e{}_minus'.format(j)])
                es_plus.append(props.loc[i, 'e{}_plus'.format(j)])
                es_mode.append(props.loc[i, 'e{}_mode'.format(j)])
                es_68.append(props.loc[i, 'e{}_68'.format(j)])

                insols.append(props.loc[i, 'insol{}'.format(j)])
                insols_med.append(props.loc[i, 'insol{}_med'.format(j)])
                insols_minus.append(props.loc[i, 'insol{}_minus'.format(j)])
                insols_plus.append(props.loc[i, 'insol{}_plus'.format(j)])

                teqs.append(props.loc[i, 'teq{}'.format(j)])
                teqs_med.append(props.loc[i, 'teq{}_med'.format(j)])
                teqs_minus.append(props.loc[i, 'teq{}_minus'.format(j)])
                teqs_plus.append(props.loc[i, 'teq{}_plus'.format(j)])

                post_paths.append('/data/user/lrosenth/legacy/final_run/{}/{}_pchains.csv'.format(starname, starname))

    # Record the planet status of each individual signal, from system_props.
    planet_status = ['' for x in range(len(masses))]
    num_planets = props['num_planets']
    index = 0
    for i in np.arange(num_stars):
        for n in np.arange(num_planets[i]):
            planet_status[index] = props.loc[i, 'status{}'.format(n+1)]
            index += 1

    # Record the posterior index of each signal.
    pl_index = np.ones(len(masses))
    index = 0
    for i in np.arange(num_stars):
        for n in np.arange(num_planets[i]):
            pl_index[index] = n + 1#props.loc[i, 'pl_index{}'.format(n+1)]
            index += 1

    planets_dict = {'hostname':starnames, 'status':planet_status, 'pl_index':pl_index,
                    'mass':masses, 'mass_med':masses_med, 'mass_minus':masses_minus,
                    'mass_plus':masses_plus, 'axis':axes, 'axis_med':axes_med,
                    'axis_minus':axes_minus, 'axis_plus':axes_plus, 'per':pers,
                    'per_med':pers_med, 'per_minus':pers_minus, 'per_plus':pers_plus,
                    'tc':tcs, 'tc_med':tcs_med, 'tc_minus':tcs_minus, 'tc_plus':tcs_plus,
                    'w':ws, 'w_med':ws_med, 'w_minus':ws_minus, 'w_plus':ws_plus,
                    'k':ks, 'k_med':ks_med, 'k_minus':ks_minus, 'k_plus':ks_plus,
                    'e':es, 'e_med':es_med, 'e_minus':es_minus, 'e_plus':es_plus,
                    'e_mode':es_mode, 'e_68':es_68, 'insol':insols,
                    'insol_med':insols_med, 'insols_minus':insols_minus,
                    'insols_plus':insols_plus, 'teq':teqs, 'teq_med':teqs_med,
                    'teq_minus':teqs_minus, 'teq_plus':teqs_plus,
                    'post_path': post_paths}
    planets_db = pd.DataFrame(planets_dict)
    planets_db.to_csv('planet_list.csv')


def plot_new_multis(planet_list_filename):
    planets = pd.read_csv(planet_list_filename)

    sys_names = planets.hostname.unique()
    new_multi_names = []
    new_multi_db = pd.DataFrame(columns=planets.columns)
    for name in sys_names:
        sys = planets.loc[planets.hostname == name].reset_index()
        sys_filtered = sys.loc[sys.status != 'N'].reset_index()
        if 'C' in sys.status.unique() and len(sys_filtered) > 1:
            new_multi_names.append(name)
            new_multi_db = new_multi_db.append(sys_filtered, ignore_index=True)

    fig, ax = plt.subplots()
    plt.title('Candidate multis')
    for name in sys_names:
        multi = new_multi_db.loc[new_multi_db.hostname == name]
        multi = multi.sort_values(by='a')
        #pdb.set_trace()
        #for i in np.arange(len(multi)):
        #    if multi.loc[i, 'status'] == 'K':
        #        ax.scatter(a, M, alpha=0.75, c='b')
        #    else:
        #        ax.scatter(a, M, alpha=0.75, c='g')
        a = multi['a']
        M = multi['mass']
        ax.plot(a, M, alpha=0.25)
        ax.scatter(a, M, alpha=0.75)#, label='Known planets')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Semi-major axis (AU)')
    ax.set_ylabel(r'Msini ($M_{Jup}$)')
    fig.savefig('mass_axis_new_multis.pdf')

    fig, ax = plt.subplots()
    plt.title('Candidate multis')
    for name in sys_names:
        multi = new_multi_db.loc[new_multi_db.hostname == name]
        multi = multi.sort_values(by='a')
        #for i in np.arange(len(multi)):
        #    if multi.loc[i, 'status'] == 'K':
        #        ax.scatter(per, M, alpha=0.75, c='b')
        #    else:
        #        ax.scatter(per, M, alpha=0.75, c='g')
        per = multi['per']
        M = multi['mass']
        ax.plot(per, M, alpha=0.25)
        ax.scatter(per, M, alpha=0.75)#, label='Known planets')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Period (days)')
    ax.set_ylabel(r'Msini ($M_{Jup}$)')
    fig.savefig('mass_period_new_multis.pdf')
