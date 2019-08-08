# Make a scatter plot of planet mass and semi-major axis.
import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import radvel

def m_a_filter(props_filename, save=True):
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
    es           = []
    es_med       = []
    es_minus     = []
    es_plus      = []

    for i in np.arange(num_stars):
        num_planets = props.loc[i, 'num_planets']
        if num_planets != 0:
            for j in np.arange(1, num_planets+1):
                starnames.append(props.loc[i, 'name'])

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

                es.append(props.loc[i, 'e{}'.format(j)])
                es_med.append(props.loc[i, 'e{}_med'.format(j)])
                es_minus.append(props.loc[i, 'e{}_minus'.format(j)])
                es_plus.append(props.loc[i, 'e{}_plus'.format(j)])
                #e = props.loc[i, 'secosw{}'.format(j)]**2 + \
                #    props.loc[i, 'sesinw{}'.format(j)]**2

    # If NASA confirmed planet numbers are in props, list
    planet_status = ['' for x in range(len(masses))]
    num_planets = props['num_planets']
    index = 0
    for i in np.arange(num_stars):
        # Identify last-found planets. Probably the new ones.
        for n in np.arange(num_planets[i]):
            planet_status[index] = props.loc[i, 'status{}'.format(n+1)]
            index += 1

    new_M = []
    old_M = []
    str_M = []
    act_M = []
    bad_M = []
    new_a = []
    old_a = []
    str_a = []
    act_a = []
    bad_a = []
    new_p = []
    old_p = []
    bad_p = []
    new_e = []
    old_e = []
    bad_e = []

    masses_k       = []
    masses_med_k   = []
    masses_minus_k = []
    masses_plus_k  = []
    axes_k         = []
    axes_med_k     = []
    axes_minus_k   = []
    axes_plus_k    = []
    pers_k         = []
    pers_med_k     = []
    pers_minus_k   = []
    pers_plus_k    = []
    es_k           = []
    es_med_k       = []
    es_minus_k     = []
    es_plus_k      = []

    for i in np.arange(len(planet_status)):
        if planet_status[i] == 'K':
            old_M.append(masses[i])
            old_a.append(axes[i])
            old_p.append(pers[i])
            old_e.append(es[i])
        elif planet_status[i] == 'C':
            new_M.append(masses[i])
            new_a.append(axes[i])
            new_p.append(pers[i])
            new_e.append(es[i])
        else:
            bad_M.append(masses[i])
            bad_a.append(axes[i])
            bad_p.append(pers[i])
            bad_e.append(es[i])

    '''
    fig, ax = plt.subplots()
    plt.title('Mass and semi-major axis scatter')
    #ax.scatter(axes, masses, c=colors, alpha=0.75)
    ax.scatter(old_a, old_M, c='b', alpha=0.75, label='Known planets')
    ax.scatter(new_a, new_M, c='g', alpha=0.75, label='Candidates')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Semi-major axis (AU)')
    ax.set_ylabel(r'Msini ($M_{Jup}$)')
    ax.legend()
    if save:
        fig.savefig('mass_axis_filter.pdf')

    fig, ax = plt.subplots()
    plt.title('Mass and period scatter')
    ax.scatter(old_p, old_M, c='b', alpha=0.75, label='Known planets')
    ax.scatter(new_p, new_M, c='g', alpha=0.75, label='Candidates')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Period (days)')
    ax.set_ylabel(r'Msini ($M_{Jup}$)')
    ax.legend()
    if save:
        fig.savefig('mass_period_filter.pdf')
    '''
    planets_dict = {'hostname':starnames, 'status':planet_status,
                    'mass':masses, 'mass_med':masses_med, 'mass_minus':masses_minus, 'mass_plus':masses_plus,
                    'axis':axes, 'axis_med':axes_med, 'axis_minus':axes_minus, 'axis_plus':axes_plus,
                    'per':pers, 'per_med':pers_med, 'per_minus':pers_minus, 'per_plus':pers_plus,
                    'k':ks, 'k_med':ks_med, 'k_minus':ks_minus, 'k_plus':ks_plus,
                    'e':es, 'e_med':es_med, 'e_minus':es_minus, 'e_plus':es_plus}
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
