import pdb

import numpy as np
import pandas as pd
import radvel

def scrape(starlist, star_db_name=None, filename='system_props.csv', fancy=True):
    """Take data from completed searches and compile into one dataframe.

    Args:
        starlist (list): List of starnames to access in current directory
        star_db_name (string [optional]): Filename of star properties dataframe
        filename (string): Path to which to save dataframe

    Note:
        If specified, compute planet masses and semi-major axes.

    """
    all_params = []
    nplanets = []

    for star in starlist:
        print(star)
        params = dict()
        params['name'] = star
        try:
            post = radvel.posterior.load(star+'/post_final.pkl')
        except (RuntimeError, FileNotFoundError):
            print('Not done looking for planets around {} yet, \
                                try again later.'.format(star))
            continue
        if fancy:
            try:
                chains = pd.read_csv(star+'/chains.csv.tar.bz2')
                masschain = np.random.normal(0, 1, len(chains))
            except (RuntimeError, FileNotFoundError):
                chains = 'empty'

        if post.params.num_planets == 1:
            if post.params['k1'].value == 0.:
                num_planets = 0
            else:
                num_planets = 1
            nplanets.append(num_planets)
        else:
            num_planets = post.params.num_planets
            nplanets.append(num_planets)
        params['num_planets'] = num_planets

        for k in post.params.keys():
            params[k] = post.params[k].value
            if fancy:
                if isinstance(chains, pd.DataFrame):
                    #pdb.set_trace()
                    params[k+'_med']   = np.median(chains[k])
                    params[k+'_minus'] = np.percentile(chains[k], 15.9)
                    params[k+'_plus']  = np.percentile(chains[k], 84.1)
        all_params.append(params)

    # Save radvel parameters as a pandas dataframe.
    props = pd.DataFrame(all_params)

    if star_db_name is not None:
        print(star)
        try:
            star_db = pd.read_csv(star_db_name)
        except (RuntimeError, FileNotFoundError):
            print('That is not a readable table. Try again.')

        # Add enough columns to account for system with the most signals.
        max_num_planets = np.amax(nplanets)
        for n in np.arange(1, max_num_planets+1):
            props['Mstar'] = np.nan
            props['M{}'.format(n)] = np.nan
            props['a{}'.format(n)] = np.nan

        # Save median star mass, uncertainties
        for star in starlist:
            try:
                props_index = props.index[props['name'] == str(star)][0]
                star_index = star_db.index[star_db['name'] == str(star)][0]
            except IndexError:
                continue
            # Save stellar mass, to be used in mass and orbital calculations.
            Mstar = star_db.loc[star_index, 'iso_mass']
            props.loc[props_index, 'Mstar'] = Mstar
            # Save stellar mass uncertainty.
            uMstar = np.mean([star_db.loc[star_index, 'iso_mass_err1'],
                              star_db.loc[star_index, 'iso_mass_err2']])
            props.loc[props_index, 'uMstar'] = uMstar

            #Make a fake posterior for stellar mass.
            if fancy:
                try:
                    chains = pd.read_csv(star+'/chains.csv.tar.bz2')
                    masschain = np.random.normal(Mstar, uMstar, len(chains))
                except (RuntimeError, FileNotFoundError):
                    chains = 'empty'
                #masschain = Mstar

            # For each found planet, compute mass and semi-major axis
            if props.loc[props_index, 'num_planets'] != 0:
                for n in np.arange(1, props.loc[props_index, 'num_planets']+1):
                    K = props.loc[props_index, 'k{}'.format(n)]
                    P = props.loc[props_index, 'per{}'.format(n)]
                    e = props.loc[props_index, 'secosw{}'.format(n)]**2 + \
                        props.loc[props_index, 'sesinw{}'.format(n)]**2
                    props.loc[props_index, 'M{}'.format(n)] = \
                        radvel.utils.Msini(K, P, Mstar, e, Msini_units='jupiter')
                    props.loc[props_index, 'a{}'.format(n)] = \
                        radvel.utils.semi_major_axis(P, Mstar)
                    if fancy:
                        if isinstance(chains, pd.DataFrame):
                            Mchain = radvel.utils.Msini(chains['k{}'.format(n)],
                                chains['per{}'.format(n)], masschain,
                                chains['e{}'.format(n)], Msini_units='jupiter')
                            props.loc[props_index, 'M{}_med'.format(n)] = \
                                np.median(Mchain)
                            props.loc[props_index, 'M{}_minus'.format(n)] = \
                                np.percentile(Mchain, 15.9)
                            props.loc[props_index, 'M{}_plus'.format(n)] = \
                                np.percentile(Mchain, 84.1)

                            achain = radvel.utils.semi_major_axis(chains['per{}'.format(n)],
                                                                  masschain)
                            props.loc[props_index, 'a{}_med'.format(n)] = \
                                np.median(achain)
                            props.loc[props_index, 'a{}_minus'.format(n)] = \
                                np.percentile(achain, 15.9)
                            props.loc[props_index, 'a{}_plus'.format(n)] = \
                                np.percentile(achain, 84.1)

            props.to_csv('system_props.csv')
    return props
