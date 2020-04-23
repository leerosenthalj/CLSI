import pdb

import numpy as np
import pandas as pd
import radvel

def insolate(T, R, a):
    return (T/5778)**4 * R**2 * a**-2

def tequil(S, alb=0.3):
    return S**-0.25 * ((1-alb)/4.)**0.25

def scrape(starlist, star_db_name=None, filename='system_props.csv', fancy=True):
    """Take data from completed searches and compile into one dataframe.

    Args:
        starlist (list): List of starnames to access in current directory
        star_db_name (string [optional]): Filename of star properties dataframe
        filename (string): Path to which to save dataframe

    Note:
        If specified, compute planet masses and semi-major axes.

    """
    nplanets   = []
    all_params = []
    all_stats  = []

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
                # Read in only last 10,000 steps.
                chains = pd.read_csv(star+'/chains.csv.tar.bz2')
                if len(chains) > 5*10**4:
                    chains = chains[-50000:]
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

        # Save quantiles and modes for all fitting-basis parameters, and for e.
        for k in post.params.keys():
            if post.params[k].vary:
                params[k] = post.params[k].value
                if fancy:
                    if isinstance(chains, pd.DataFrame):
                        params[k+'_med']   = np.median(chains[k])
                        params[k+'_minus'] = np.percentile(chains[k], 15.9)
                        params[k+'_plus']  = np.percentile(chains[k], 84.1)
                        hist, bin_edges    = np.histogram(chains[k], bins=200)
                        params[k+'_mode']  = bin_edges[np.argmax(hist)]
        if num_planets > 0:
            for n in np.arange(1, num_planets+1):
                ekey = 'e{}'.format(n)
                wkey = 'w{}'.format(n)
                tkey = 'tp{}'.format(n)
                params[ekey] = post.params[ekey].value
                params[wkey] = post.params[wkey].value
                params[tkey] = post.params[tkey].value
                if fancy:
                    if isinstance(chains, pd.DataFrame):
                        #echain = chains['secosw{}'.format(n)]**2 + \
                        #         chains['sesinw{}'.format(n)]**2
                        echain = chains[ekey]
                        params[ekey+'_med']   = np.median(echain)
                        params[ekey+'_minus'] = np.percentile(echain, 15.9)
                        params[ekey+'_plus']  = np.percentile(echain, 84.1)
                        hist, bin_edges       = np.histogram(echain, bins=200,
                                                             range=(0, 1))
                        params[ekey+'_mode'] = bin_edges[np.argmax(hist)]
                        params[ekey+'_68'] = np.percentile(echain, 68.2)
                        std = 0.5*(params[ekey+'_plus'] - params[ekey+'_minus'])
                        if params[ekey+'_mode'] > 2*std:
                            params[ekey+'_stat'] = 'e'
                        else:
                            params[ekey+'_stat'] = 'c'

        all_params.append(params)

        # Collect observation stats. Nobs, baseline, median error for each inst.
        stats = dict()
        stats['name'] = star
        for like in post.likelihood.like_list:
            tel = like.telvec[0]
            stats['Nobs_'+tel] = len(like.x)
            stats['baseline_'+tel] = np.amax(like.x) - np.amin(like.x)
            stats['med_err_'+tel] = np.median(like.yerr)
        stats['Nobs'] = len(post.likelihood.x)
        stats['baseline'] = np.amax(post.likelihood.x) - np.amin(post.likelihood.x)

        all_stats.append(stats)


    # Save radvel parameters as a pandas dataframe.
    props = pd.DataFrame(all_params)
    props.to_csv('system_props_no_mass.csv')

    all_stats_db = pd.DataFrame(all_stats)
    all_stats_db.to_csv('observation_stats.csv')


    if star_db_name is not None:
        try:
            star_db = pd.read_csv(star_db_name)
        except (RuntimeError, FileNotFoundError):
            print('That is not a readable table. Try again.')

        # Add enough columns to account for system with the most signals.
        max_num_planets = np.amax(nplanets)
        props['Mstar']  = np.nan
        props['uMstar'] = np.nan
        for n in np.arange(1, max_num_planets+1):
            props['M{}'.format(n)] = np.nan
            props['a{}'.format(n)] = np.nan

        merge = props.merge(star_db, on='name')
        props = merge
        props.to_csv('system_props_specmatch.csv')
        # Save median star mass, uncertainties
        for index, row in props.iterrows():
            star = props.loc[index, 'name']
            print(star)

            # Save stellar params and errs, to be used in derived calculations.
            # IMPORTANT: STANDARDIZE AFTER MERGING SYNTH AND EMP PARAMS. 10/25
            Mstar  = props.loc[index, 'mass_c']
            uMstar = props.loc[index, 'mass_err_c']
            Rstar  = props.loc[index, 'radius_c']
            uRstar  = props.loc[index, 'radius_err_c']
            Tstar  = props.loc[index, 'teff_c']
            uTstar = props.loc[index, 'teff_err_c']

            props.loc[index, 'Mstar']  = Mstar
            props.loc[index, 'uMstar'] = uMstar

            #Make a fake posterior for stellar mass.
            if fancy:
                try:
                    chains = pd.read_csv(star+'/chains.csv.tar.bz2')
                    if len(chains) > 5*10**4:
                        chains = chains[-50000:]
                except (RuntimeError, FileNotFoundError):
                    chains = 'empty'
                if ~np.isnan(Mstar): #try:
                    masschain = np.random.normal(Mstar, uMstar, len(chains))
                    radchain  = np.random.normal(Rstar, uRstar, len(chains))
                    tempchain = np.random.normal(Tstar, uTstar, len(chains))
                else: #except (RuntimeError, ValueError):
                    masschain = np.random.normal(1., 0.05, len(chains))
                    radchain  = np.random.normal(1., 0.05, len(chains))
                    tempchain = np.random.normal(5700., 100., len(chains))
                    print('BAD')

            # If reading posteriors, make dictionary for physical param chains.
            if fancy:
                pdict = {}

            # For each found planet, compute mass and semi-major axis
            if props.loc[index, 'num_planets'] != 0:
                for n in np.arange(1, props.loc[index, 'num_planets']+1):
                    K = props.loc[index, 'k{}'.format(n)]
                    P = props.loc[index, 'per{}'.format(n)]
                    e = props.loc[index, 'secosw{}'.format(n)]**2 + \
                        props.loc[index, 'sesinw{}'.format(n)]**2

                    props.loc[index, 'M{}'.format(n)] = \
                        radvel.utils.Msini(K, P, Mstar, e, Msini_units='jupiter')
                    props.loc[index, 'a{}'.format(n)] = \
                        radvel.utils.semi_major_axis(P, Mstar)

                    if fancy:
                        if isinstance(chains, pd.DataFrame):
                            echain = chains['secosw{}'.format(n)]**2 + \
                                     chains['sesinw{}'.format(n)]**2
                            wchain = np.arctan(chains['sesinw{}'.format(n)]/
                                               chains['secosw{}'.format(n)])
                            Mchain = radvel.utils.Msini(chains['k{}'.format(n)],
                                chains['per{}'.format(n)], masschain,
                                echain, Msini_units='jupiter')
                            achain = radvel.utils.semi_major_axis(chains[
                                                  'per{}'.format(n)], masschain)
                            insolchain = insolate(tempchain, radchain, achain)
                            teqchain   = tequil(insolchain)
                            # Save physical chains.
                            # M, a, e, w, K, P, tc, tp, insol, teq
                            pdict['M{}'.format(n)]     = Mchain
                            pdict['a{}'.format(n)]     = achain
                            pdict['e{}'.format(n)]     = echain
                            pdict['w{}'.format(n)]     = wchain
                            pdict['tp{}'.format(n)]    = chains['tc{}'.format(n)]
                            pdict['k{}'.format(n)]     = chains['k{}'.format(n)]
                            pdict['per{}'.format(n)]   = chains['per{}'.format(n)]
                            pdict['tc{}'.format(n)]    = chains['tc{}'.format(n)]
                            pdict['insol{}'.format(n)] = insolchain
                            pdict['teq{}'.format(n)]   = teqchain

                            # Save fitting and physical quantiles and modes.
                            props.loc[index, 'M{}_med'.format(n)] = \
                                np.median(Mchain[~np.isnan(Mchain)])
                            props.loc[index, 'M{}_minus'.format(n)] = \
                                np.percentile(Mchain[~np.isnan(Mchain)], 15.9)
                            props.loc[index, 'M{}_plus'.format(n)] = \
                                np.percentile(Mchain[~np.isnan(Mchain)], 84.1)

                            props.loc[index, 'a{}_med'.format(n)] = \
                                np.median(achain[~np.isnan(achain)])
                            props.loc[index, 'a{}_minus'.format(n)] = \
                                np.percentile(achain[~np.isnan(achain)], 15.9)
                            props.loc[index, 'a{}_plus'.format(n)] = \
                                np.percentile(achain[~np.isnan(achain)], 84.1)

                            props.loc[index, 'w{}_med'.format(n)] = \
                                np.median(pdict['w{}'.format(n)])
                            props.loc[index, 'w{}_minus'.format(n)] = \
                                np.percentile(pdict['w{}'.format(n)], 15.9)
                            props.loc[index, 'w{}_plus'.format(n)] = \
                                np.percentile(pdict['w{}'.format(n)], 84.1)

                            props.loc[index, 'tp{}_med'.format(n)] = \
                                np.median(pdict['tp{}'.format(n)])
                            props.loc[index, 'tp{}_minus'.format(n)] = \
                                np.percentile(pdict['tp{}'.format(n)], 15.9)
                            props.loc[index, 'tp{}_plus'.format(n)] = \
                                np.percentile(pdict['tp{}'.format(n)], 84.1)

                            props.loc[index, 'insol{}_med'.format(n)] = \
                                np.median(insolchain[~np.isnan(insolchain)])
                            props.loc[index, 'insol{}'.format(n)] = \
                                insolate(Tstar, Rstar, props.loc[index, 'a{}'.format(n)])
                            props.loc[index, 'insol{}_minus'.format(n)] = \
                                np.percentile(insolchain[~np.isnan(insolchain)], 15.9)
                            props.loc[index, 'insol{}_plus'.format(n)] = \
                                np.percentile(insolchain[~np.isnan(insolchain)], 84.1)

                            props.loc[index, 'teq{}_med'.format(n)] = \
                                np.median(teqchain[~np.isnan(teqchain)])
                            props.loc[index, 'teq{}'.format(n)] = \
                                tequil(props.loc[index, 'insol{}'.format(n)])
                            props.loc[index, 'teq{}_minus'.format(n)] = \
                                np.percentile(insolchain[~np.isnan(insolchain)], 15.9)
                            props.loc[index, 'teq{}_plus'.format(n)] = \
                                np.percentile(insolchain[~np.isnan(insolchain)], 84.1)

                # Save star's physical, thinned chain.
                if fancy:
                    pchains = pd.DataFrame.from_dict(pdict)
                    # Get rid of any rows with nans due to negative mstar sample.
                    pchains = pchains.loc[~np.isnan(Mchain)]
                    # Thin the chains. Not doing for now, since only saving final steps.
                    #pchains = pchains.iloc[::10, :]
                    pchains.to_csv(star+'/{}_pchains.csv'.format(star))

            props.to_csv('system_props.csv')
    return props
