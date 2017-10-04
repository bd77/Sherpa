# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:15:49 2017
Script used to produce the figures of the Atlas, starting from
the excel resutls file from Bart,
and the nc files for Bart.
@author: peduzem
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import numpy as np

import pandas as pd

from sherpa_globals import path_model_cdf_test


def figsize(scale):
    '''Square figure size'''
    fig_width_pt = 390  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
#    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
    fig_height = fig_width #* golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size

def figsizer(scale): # rectangular figure
    '''Rectangualar figure size'''
    fig_width_pt = 390  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
    fig_height = fig_width * golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size

if __name__ == '__main__':

    # -------------------------------------------------------------------------
    # Inputs:
    # results path
    path_results = 'D:\\sherpa.git\\Sherpa\\atlas\\'
    # path for Figure outputs
    path_figures = path_results + 'figures\\'
    # results filename from Bart
    filename = 'FUA112_source_apportionment_20170904.xlsx'

    # ASCI names and URAU code fo cities (from Marco)
    fuanames = pd.read_excel((path_results+'fua_asci_uraucodes2.xlsx'),
                             index_col=1)

    # Other inputs:
    # Cities with which are big enough to be displayed with cores:
    wcore = pd.read_excel(path_results+'onlyfua_city-fua.xls',
                          sheetname='withcore', header=[0], index_col=[0],
                          na_values='nan').sortlevel()
    # precursor list
    precursor_lst = ['NOx', 'NMVOC', 'NH3', 'PPM', 'SOx']
    # sector aggregation:
    aggr_dict = {'Industry': [1, 3, 4], 'Residential': [2],
                 'Agriculture': [10], 'Transport': [7], 'Other': [5, 6, 8, 9]}
    sect_aggr = list(aggr_dict.keys())

    # -------------------------------------------------------------------------
    # Get data from Bart's results
    df = pd.read_excel((
            path_results + filename),
            sheetname='chimere', skiprows=3, header=[1],
            index_col=[2, 3, 5, 4], na_values='nan').sortlevel()
    df.sort_index(inplace=True)
    # Cities in the results file
    cities = list(set(df.index.get_level_values('target')))

    # -------------------------------------------------------------------------
    # Get the emission data and put it in a dic of multindex dataframes
    # Dictionary of dataframes:
    emis_ag = {}
    # Define multi index for each dataframe:
    it_ag = [precursor_lst, sect_aggr]
    in_ag = pd.MultiIndex.from_product(it_ag, names=['precursor', 'aggS'])

    # loop ofer all the cities
    for city in cities:
        # Cread dataframe
        emis_ag[city] = pd.DataFrame(index=in_ag,
                                     columns=['core', 'comm', 'fua'])
        # Sort df, it is important for copying values in the right place
        emis_ag[city].sort_index(inplace=True)
        # Fill dataframe with aggregated values
        for key in aggr_dict.keys():
            emis_ag[city].loc[(precursor_lst, key), 'core'] = (
                    df['DE'].loc[city, (city+'_City'), :, aggr_dict[key]].sum(
                            skipna=True, level=2).values*2)
            # not all city have a commuting zone, if there is no commuting zone
            # emission value is zero
            try:
                emis_ag[city].loc[(precursor_lst, key), 'comm'] = (
                    df['DE'].loc[city, (city+'_Comm'), :, aggr_dict[key]].sum(
                             skipna=True, level=2).values*2)
            except(KeyError):  # for places without a commuting zone
                emis_ag[city].loc[(precursor_lst, key), 'comm'] = 0
        # Emissions of Fua are the sum of city and commuting zone
        emis_ag[city]['fua'] = emis_ag[city].sum(axis=1)

    # -------------------------------------------------------------------------
    # Create polar polar plots for all precursor and all cities:
    # Prepare data
    index = aggr_dict.keys()
    plt.close('all')
    for prec in precursor_lst:
        for city in cities:
            print(city)
#            city = 'Berlin'
            # list of dataframe from which to make plots
            df_polar = []
            # liverpool does not have a commuting zone, only the center
            # should be considered
            if city is not 'Liverpool':
                # get data for FUA
                df_fua = (emis_ag[city]['fua'].loc[
                          prec, sect_aggr].unstack(level=1).T / 1000)
                df_fua = df_fua.reindex(sect_aggr)
                df_polar.append([df_fua.fillna(value=0)])
                lableg = ['Greater city', 'City']
                colors = ['blue', 'red']
                colorsf = ['None', 'red']
                hatch = ['', 3*'///']
            if city is 'Liverpool':
                lableg = ['Greater city']
                colors = ['red']
                colorsf = ['red']
                hatch = [3*'///']
            # if city is big enough to consider the core separately
            if city in wcore.index:
                df_city = (emis_ag[city]['core'].loc[
                           prec, sect_aggr].unstack(level=1).T / 1000)
                df_city = df_city.reindex(sect_aggr)
                df_polar.append([df_city.fillna(value=0)])
    # -------------------------------------------------------------------------
    # Create polar plots
            # fontsize:
            fts = 5
            # create figure
            fig = plt.figure(figsize=figsize(0.3))
            ax = fig.add_subplot(111, projection="polar")
            # create grid
            ax.grid(True)
            ax.yaxis.grid(color='#aab0b7', lw=0.5)
            ax.xaxis.grid(color='#aab0b7', lw=0.5)
            # set bored of figure area
            for spine in ax.spines.values():
                spine.set_edgecolor('#aab0b7')
                spine.set_zorder(0)
                spine.set_linewidth(0.5)

            # Define angles of y axisi in polar plots
            theta = (np.arange(len(df_polar[0][0])) /
                     float(len(df_polar[0][0]))*2.*np.pi)

            # Define tick labels for y axis
            ax.set_xticks(theta)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_horizontalalignment('center')
                tick.label1.set_verticalalignment('top')
            # get radial labels away from plotted line (set angle)
            ax.set_rlabel_position(90)

            # plot data
            plots = []
            mpl.rcParams['hatch.color'] = 'red'
            for it, dfdata in enumerate(df_polar):
                line, = ax.plot(theta, df_polar[it][0], color=colors[it],
                                label=None, zorder=3, lw=0.7)
                ax.fill(theta, df_polar[it][0], color=colorsf[it], alpha=0.3,
                        zorder=3, lw=0.3)
                ax.fill(theta, df_polar[it][0], edgecolor=colorsf[it],
                        alpha=0.3, hatch=hatch[it], zorder=3, color='None',
                        lw=0.4)
                plots.append(line,)

            def _closeline(line):
                x, y = line.get_data()
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)
            [_closeline(l) for l in plots]

            # maximums y value... so I can put the label
            maxy = max(df_polar[0][0][prec])

            # set yticks so that I have 5 y ticks for each figure
            ylim_lst = [1, 2, 5, 10, 15, 20, 40, 50, 60, 80, 100]
            for lim in ylim_lst:
                if maxy <= lim:
                    ylim = lim
                    ax.set_ylim(top=ylim)
                    plt.yticks(np.arange(0, (ylim+ylim/5), ylim/5), zorder=6)
                    break
                else:
                    ylim = maxy
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)

            # Substituting the axis label with a txt label otherwise
            # the lable is below the grid (probably a bug)
            for it in np.arange(len(theta)):
                ax.text(theta[it], ylim*1.15, sect_aggr[it], va='center',
                        ha='center', fontsize=(fts+2))
            ax.set_xticklabels('')

            plt.legend(handles=plots, labels=lableg,
                       bbox_to_anchor=(1.45, 0.2), fontsize=fts, frameon=False)

            ax.tick_params(axis='y', labelsize=fts, zorder=4)

            plt.title("{} \n [kton/y]".format(prec), fontsize=(fts+2),
                      y=0.9, x=0, weight='bold')
            citynospace = city.translate({ord(c): None for c in ' \n\t\r'})

            plt.savefig((path_figures + fuanames.loc[city, 'URAU_CODE'] +
                         '_{}'.format(citynospace) +
                         '_emi_FUA_{}'.format(prec)+'.png'),
                        dpi=1000, bbox_inches='tight', transparent=True)
            plt.show()

    # -------------------------------------------------------------------------
    # Get the source allocaltion data and put it in a dic of multindex
    # dataframes
    # Dictionary of dataframes:
    dc_ag = {}  # Dictionary of dataframes to store data
    # name of states corresponding to each city
    df_names = pd.read_excel(path_results + filename,
                             sheetname='cities', index_col=[0])
    # indicator for source allocation (column name)
    ind = 'relative_potential'
    for city in cities:
#        city = 'Liverpool'
        print(city)
        # create Dataframe
        columns = ['City', 'Comm', 'National', 'International']
        dc_ag[city] = pd.DataFrame(index=sect_aggr,
                                   columns=columns)
        dc_ag[city].sort_index(inplace=True)

        # Get name of country to rename index to read for international contr.:
        for col in columns[0:4]:
#            col = 'City'
            if col == 'International':
                name = df_names['country'].loc[city]
            else:
                name = city
            for key in aggr_dict.keys():
                try:
                    dc_ag[city].loc[key, col] = sum(df[ind].loc[city,'{}_{}'.format(name, col), :, aggr_dict[key]].values)
                except(KeyError):
                    dc_ag[city][col] = 0

        dc_ag[city]['Total'] = dc_ag[city].sum(axis=1)
        if city not in wcore.index:
            dc_ag[city]['Greater city'] = dc_ag[city][['City',
                                                       'Comm']].sum(axis=1)
            del dc_ag[city]['City']
            del dc_ag[city]['Comm']
            nat = dc_ag[city]['Greater city'].sum()
            intern = nat + dc_ag[city]['National'].sum()
            dc_ag[city] = dc_ag[city][['Greater city', 'National',
                                      'International', 'Total']]
            dc_ag[city].loc['bottom'] = [0, nat, intern, 0]
            dc_ag[city].loc['no control'] = [0, 0, 0,
                                             (100-dc_ag[city]['Total'].sum())]
        else:
            ci = dc_ag[city]['City'].sum()
            nat = ci+dc_ag[city]['Comm'].sum()
            intern = nat + dc_ag[city]['National'].sum()
            dc_ag[city].loc['bottom'] = [0, ci, nat, intern, 0]
            dc_ag[city].loc['no control'] = [0, 0, 0, 0,
                                             (100-dc_ag[city]['Total'].sum())]

        plt.close('all')
        fts = 10
        f, ax1 = plt.subplots(1, figsize=figsizer(1))
        ax1.yaxis.grid(color='black')

        colors = {'PPM': 'blue', 'SOx': 'gold', 'NOx': 'red', 'NH3': 'green',
                  'NMVOC': '#black', 'Transport': 'red', 'Energy': 'blue',
                  'Industry': 'gold', 'Production': '#8B4789', 'Waste': '#00FFFF',
                  'Agriculture': 'green', 'Residential': 'blue',
                  'Offroad': '#8470ff', 'Extraction': '#00FF00',
                  'Other': 'skyblue', 'Salt': '#ccffe5', 'Dust': '#ffffcc',
                  'no control': '#dddddd', 'bottom': 'None'}

        index = []
        index.append('bottom')
        index.extend(sect_aggr)
        ncontrol = ['no control']
        index.extend(ncontrol)

        type_colors = [colors[k] for k in index]
        if city == 'Liverpool':
            del dc_ag[city]['Comm']
        dc_ag[city].reindex(index).T.plot(kind='barh', stacked=True,
                               color=type_colors, xticks=None, legend=False,
                               ax=ax1,fontsize=fts)

        handles, labels = ax1.get_legend_handles_labels()
        labels = labels[1:7]
        handles = handles[1:7]
        ax1.legend(handles, labels, fontsize=fts-4, frameon=False)
        ax1.set_xlabel('Percentage of total mass', fontsize=fts)
        ytlab = [label.get_text() for label in ax1.get_yticklabels()]
        dict_ytlab = {'Comm': 'Commuting \n Zone',
                      'National': 'Rest of \n the country'}
        for i, ytick in enumerate(ytlab):
            if ytick in dict_ytlab.keys():
                ytlab[i] = dict_ytlab[ytick]

        ax1.set_yticklabels(ytlab)
        citynospace = city.translate({ord(c): None for c in ' \n\t\r'})
        plt.savefig((path_figures+fuanames.loc[city, 'URAU_CODE'] +
                    '_{}'.format(citynospace)+'_conc_FUA.png'),
                    dpi=1000, bbox_inches='tight', transparent=True)
        plt.show()

    # -------------------------------------------------------------------------
    # Create netcdf files for bart
    from sherpa_auxiliaries_epe import gridint_toarray, write_nc
    import os as os

    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()

    # Make nc file for Bart to use
    path_results = 'D:\\sherpa.git\\Sherpa\\atlas\\ncfiles\\'
    # ASCI names and URAU code fo cities (from Marco)
    gcity_map = pd.read_table('atlas\\151_fuas_cities_urau.txt',
                              index_col=['FUA_CODE', 'ASCII_NAME'])
    fuanames = pd.read_excel(('atlas\\fua_asci_uraucodes151.xls'),
                             index_col=1)
    for city in fuanames.index:
        print('making nc file for', city)
#        city = 'Paris'
#        city2 = 'Antwerpen'
        urau_code = fuanames.loc[city]['URAU_CODE']  # URAU_CODE == FUA_CODE
#        urau_code2 = fuanames.loc[city2]['URAU_CODE']

        cityc_codes = gcity_map.loc[urau_code]['URAU_COD_2']
#        cityc_codes2= gcity_map.loc[urau_code2]['URAU_COD_2']
#        cityc_codes2= np.array(gcity_map.loc[urau_code2]['URAU_COD_2'])
        # City area:
        # create a matrix of the same size as sherpas grid
        area_core = np.zeros((1, len(lat_array), len(lon_array)))
        path_area = path_results+'{}_City.nc'.format(city)
        if not os.path.exists(path_area):
            for code in cityc_codes:
                print(code)
                area_city = gridint_toarray('GCITY_CODE', 'parea', code)
                area_core = area_core + area_city
            write_nc(area_core, path_area, 'AREA', '%', addnutsid=True)

        # Commuting zone area:
        path_area = path_results+'{}_Comm.nc'.format(city)
        if not os.path.exists(path_area):
            area_fua = gridint_toarray('FUA_CODE', 'parea', urau_code)
            area_comm = area_fua - area_core
            write_nc(area_comm, path_area, 'AREA', '%', addnutsid=True)

        # Rest of the country area
        path_area = path_results+'{}_National.nc'.format(city)
        if not os.path.exists(path_area):
            area_country = gridint_toarray('NUTS_Lv0', 'parea',urau_code[:2])
            area_national = area_country - area_fua
            write_nc(area_national, path_area, 'AREA', '%', addnutsid=True)

    # -------------------------------------------------------------------------



#    emis, emis_ag_denise = precompute(path_area)
    # BAR PLOT
#    fig = plt.figure(figsize=figsize(1))
#    ax1 = fig.add_subplot(111, facecolor='None')
#    for spine in ax1.spines.values():
#        spine.set_edgecolor('black')
#        spine.set_linewidth(0.5)
#    ax1.yaxis.grid(color='#dddddd')
#    type_colors = ['mediumpurple',  # industry
#                   'tomato',  # residential
#                   'mediumseagreen',  # transport
#                   'gold',  # agriculture
#                   '#6890F0',  # other
#                   ]
#    df1[sect_display].plot(kind='bar', stacked=False, alpha=1,
#                           color=type_colors, xticks=None, legend=False,
#                           hatch='///', ax=ax1)
#    df2[sect_display].plot(kind='bar', stacked=False, alpha=1, xticks=None,
#                           ax=ax1, color=type_colors)
#    ax1.get_xaxis().set_ticks([])
#    ax1.set_ylabel('PPM emissions [kton/year]')
#    ax1.set_xlabel('')
#    hatch_patch = mpatches.Patch(hatch='///', color='None')
#    handles, labels = ax1.get_legend_handles_labels()
#    labels = labels[5:10] + ['com. zone']
#    handles = handles[5:10]
#    handles.append(hatch_patch)
#    ax1.legend(handles, labels, facecolor='w', edgecolor='black')
#    plt.savefig('atlas//Paris.pdf')


#    ## DENISE:
#    cityname = 'PARIS'
#    level ='GCITY_CODE'
#    code = 'FR001K1'
#    path_area = 'workdir\\Paris_city12.nc'#.format(cityname,level,code)
#    if not os.path.exists(path_area):
#        areacountry = gridint_toarray(level, 'parea', code)
#        write_nc(areacountry, path_area, 'AREA', '%', addnutsid=True)
#    emis, emis_ag_denise = precompute(path_area)
#
#    level = 'FUA_CODE'
#    code = 'FR001L1'
#    path_area = 'workdir\\Paris_FUA.nc'#.format(cityname, level, code)
#    if not os.path.exists(path_area):
#        areacountry2 = gridint_toarray(level, 'parea', code)
#        write_nc(areacountry2, path_area, 'AREA', '%', addnutsid=True)
#    emis, emis_ag_denise_fua = precompute(path_area)
#    emis_ag_denise_cz = emis_ag_denise_fua - emis_ag_denise
#    path_area = 'workdir\\Paris_comm.nc' #.format(cityname, level, code)
#    write_nc((areacountry2-areacountry), path_area, 'AREA', '%', addnutsid=True)


#def plot_spider(emi_sum, prec):  # Added by Ema for spyder plot
#    """Spyder plot for emissions, only if the point is in GCITY
#    @author: peduzem
#
#    """
##    prec = 'PPM'
#    if 'GCITY_CODE' in emi_sum.columns:
#        sect_display = emi_sum['GCITY_CODE'].loc[prec].index.get_level_values(0)
#
#        index = [sect_display[1], sect_display[3], sect_display[4],
#                 sect_display[0], sect_display[2]]
#        df_city = emi_sum['GCITY_CODE'].loc[prec].reindex(index)/1000
#        fig = plt.figure()
##        fig = plt.figure(figsize=figsize(1))
#        ax = fig.add_subplot(111, projection="polar")
#        ax.grid(True)
#        ax.yaxis.grid(color='#dddddd')
#        ax.xaxis.grid(color='#dddddd')
#        colors = ['gold', 'tomato']
#        for spine in ax.spines.values():
#            spine.set_edgecolor('None')
#        theta = np.arange(len(df_city))/float(len(df_city))*2.*np.pi
#        ax.set_xticks(theta)
#        for tick in ax.yaxis.get_major_ticks():
#            tick.label1.set_horizontalalignment('center')
#            tick.label1.set_verticalalignment('top')
#
#        df = []
#
#        if 'FUA_CODE' in emi_sum.columns:
#            df_fua = emi_sum['FUA_CODE'].loc[prec].reindex(index)/1000
#            df.append([df_fua])
#        df.append([df_city])
#        plots = []
#        for it, dfdata in enumerate(df):
#            line, = ax.plot(theta, df[it][0], color=colors[it],
#                            marker="o", label=None, zorder=3)  # , zorder = -3)
#            ax.fill(theta, df[it][0], color=colors[it], alpha=0.3,
#                    zorder=3)
#            plots.append(line,)
#
#        def _closeline(line):
#            x, y = line.get_data()
#            x = np.concatenate((x, [x[0]]))
#            y = np.concatenate((y, [y[0]]))
#            line.set_data(x, y)
#        [_closeline(l) for l in [plots[0], plots[1]]]
#
#        ax.set_rlabel_position(90)  # get radial labels away from plotted line
#
#        maxy = max(max(df_fua), max(df_city))
#        for it in np.arange(len(theta)):
#            ax.text(theta[it], maxy*1.15, index[it], va='center',
#                    ha='center', fontsize=12)
#        ax.set_xticklabels('')
#
#        plt.legend(handles=[plots[0], plots[1]], labels=['FUA', 'city core'],
#                   bbox_to_anchor=(0, 1), loc='center left')
#        ax.tick_params(axis='y', labelsize=10, zorder=4)
#        plt.title("{} emissions [kton/y]".format(prec), fontsize=14, y=1.2)
#        return fig


#def precompute(path_area):
#    """
#    Precompute function: computes given an area the correpsonding emission:
#    Input:
#        - path_area: path to the area
#    Output:
#        - emi: total emission per pollutant [Mg]
#        - emi_sec: total emissions  per pollutant per sector [Mg]
#    Files needed:
#        - nc of the area of each country in:
#            'workdir\\area_{}.nc'.format(country)
#
#    Created on Thu Apr 20 16:26:42 2017
#    @author: peduzem
#    """
#    # Total emissions by country and precursor
#    # -------------------------------------------------------------------------
#    # create a dictionary with emissions per precursor, macrosector and postion
#    rootgrp = Dataset(path_model_cdf_test, 'r')
#    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
#    rootgrp.close()
#    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)
#    # Area of each cell in the domain
#    path_surf_nc = 'input/JRC01.nc'
#    rootgrp = Dataset(path_surf_nc, mode='r')
#    surf = rootgrp.variables['surface'][:]
#    # area_units = rootgrp.variables['surface'].units
#    rootgrp.close()
#
#    # Percentage of each cell in the domain
#    rootgrp = Dataset(path_area, mode='r')
#    area = rootgrp.variables['AREA'][:]
#    # area_units = rootgrp.variables['AREA'].units
#    rootgrp.close()
#    emis = pd.DataFrame(index=precursor_lst, columns=['emi'])
#
#    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
#    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
#                 'other']
#    iterables = [precursor_lst, ms_list]
#    index = pd.MultiIndex.from_product(iterables, names=['precursor', 'MS'])
#    emis_sec = pd.DataFrame(index=index, columns=['emi'])
#
#    it_ag = [precursor_lst, sect_aggr]
#    in_ag = pd.MultiIndex.from_product(it_ag, names=['precursor', 'aggS'])
#    emis_ag = pd.DataFrame(index=in_ag, columns=['emi'])
#
#    emis['emi'] = [np.sum(np.sum((emission_dict[precursor] *
#                   surf * area / 100), axis=0)) for precursor in precursor_lst]
#    emis_sec['emi'] = [np.sum(emission_dict[precursor][ms] * surf * area / 100)
#                       for precursor in precursor_lst
#                       for ms in np.arange(0, 10)]
#    emis_ag.sort_index(inplace=True)
#    emis_sec.sort_index(inplace=True)
#
#    emis_ag['emi'].loc[:, sect_aggr[0]] = np.sum(
#           emis_sec['emi'].loc[:, ms].values for ms in ['MS1', 'MS3', 'MS4'])
#    emis_ag['emi'].loc[:,sect_aggr[1]] = np.sum(emis_sec['emi'].loc[:,ms].values for ms in ['MS2'])
#    emis_ag['emi'].loc[:,sect_aggr[2]] = np.sum(emis_sec['emi'].loc[:,ms].values for ms in ['MS10'])
#    emis_ag['emi'].loc[:,sect_aggr[3]] = np.sum(emis_sec['emi'].loc[:,ms].values for ms in ['MS7'])
#    emis_ag['emi'].loc[:,sect_aggr[4]] = np.sum(emis_sec['emi'].loc[:,ms].values for ms in ['MS5','MS6', 'MS8', 'MS9'])
#    return emis, emis_ag
#    emis_ag['emi'].loc[:,sect_aggr[0]] = [emis_sec['emi'].loc[precursor,'MS1'] for precursor in precursor_lst]
#    pd.merge(emis_ag['emi'], emis_sec['emi'], on='emi', how='outer')
#    pd.merge(emis_ag['emi'].loc[:,sect_aggr[0], emis_sec['emi'].loc[:,'MS1'], on='emi', how='outer')
#    emis_sec['emi'].loc[:,'MS1']
#    return emi, emi_sec
