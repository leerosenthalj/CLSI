# Make a scatter plot of planet mass and semi-major axis.
import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import radvel

def m_a_scatter(props_filename, save=True):
    props = pd.read_csv(props_filename)

    # Collect the starname/mass/axis sets for each planet.
    num_stars = len(props)
    starnames = []
    masses = []
    axes = []
    pers = []

    for i in np.arange(num_stars):
        num_planets = props.loc[i, 'num_planets']
        if num_planets != 0:
            for j in np.arange(1, num_planets+1):
                starnames.append(props.loc[i, 'name'])
                masses.append(props.loc[i, 'M{}'.format(j)])
                axes.append(props.loc[i, 'a{}'.format(j)])
                pers.append(props.loc[i, 'per{}'.format(j)])

    # If NASA confirmed planet numbers are in props, list
    new_planet = np.zeros(len(masses))
    index = 0
    if 'nasa_num_planets' in props:
        new = props['num_planets'] - props['nasa_num_planets']
        index = 0
        for i in np.arange(num_stars):
            # Identify last-found planets. Probably the new ones.
            if new[i] > 0:
                for n in np.arange(props.loc[i, 'num_planets']):
                    if n+1 > props.loc[i, 'nasa_num_planets']:
                        new_planet[index] = 1
                    index += 1
            elif props.loc[i, 'num_planets'] == 0:
                pass
            else:
                index += props.loc[i, 'num_planets']
        #pdb.set_trace()
    '''
    if 'nasa_num_planets' in props:
        for i in np.arange(num_stars):
            num_planets = props.loc[i, 'num_planets']
            nasa_num_planets = props.loc[i, 'nasa_num_planets']
            if num_planets != nasa_num_planets:
                for n in np.arange(1,num_planets+1):
                    new_planet[index] = 1
                    index += 1
            else:
                index += num_planets
    '''

    new_M = []
    old_M = []
    new_a = []
    old_a = []
    new_p = []
    old_p = []

    for i in np.arange(len(new_planet)):
        if new_planet[i] == 0:
            old_M.append(masses[i])
            old_a.append(axes[i])
            old_p.append(pers[i])
        else:
            new_M.append(masses[i])
            new_a.append(axes[i])
            new_p.append(pers[i])
    #pdb.set_trace()

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
        fig.savefig('mass_axis.pdf')
        #planets_dict = {'hostname':starnames, 'mass':masses, 'a':axes}
        #planets_db = pd.Dataframe(planets_dict)
        #planets_db.to_csv('planet_list.csv')

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
        fig.savefig('mass_period.pdf')


def m_a_filter(props_filename, save=True):
    props = pd.read_csv(props_filename)

    # Collect the starname/mass/axis sets for each planet.
    num_stars = len(props)
    starnames = []
    masses    = []
    axes      = []
    pers      = []
    es        = []

    for i in np.arange(num_stars):
        num_planets = props.loc[i, 'num_planets']
        verdict = props.loc[i, 'Verdict']
        if num_planets != 0:
            for j in np.arange(1, num_planets+1):
                starnames.append(props.loc[i, 'name'])
                masses.append(props.loc[i, 'M{}'.format(j)])
                axes.append(props.loc[i, 'a{}'.format(j)])
                pers.append(props.loc[i, 'per{}'.format(j)])
                e = props.loc[i, 'secosw{}'.format(j)]**2 + \
                    props.loc[i, 'sesinw{}'.format(j)]**2
                es.append(e)

    # If NASA confirmed planet numbers are in props, list
    planet_status = ['' for x in range(len(masses))]
    if 'nasa_num_planets' in props:
        num_planets = props['num_planets']
        index = 0
        for i in np.arange(num_stars):
            # Identify last-found planets. Probably the new ones.
            for n in np.arange(num_planets[i]):
                planet_status[index] = props.loc[i, 'status{}'.format(n+1)]
                index += 1

    new_M = []
    old_M = []
    bad_M = []
    new_a = []
    old_a = []
    bad_a = []
    new_p = []
    old_p = []
    bad_p = []
    new_e = []
    old_e = []
    bad_e = []

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

        planets_dict = {'hostname':starnames, 'mass':masses, 'a':axes,
                        'per':pers, 'e':es, 'status':planet_status}
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
