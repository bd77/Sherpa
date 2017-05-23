# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:26:42 2016
@author: peduzem

Script to calculate the blame matrix

"""
# for importing matlab files
import scipy.io as sio
from PIL import Image
from osgeo import gdal, ogr, osr  #conda install -c conda-forge gdal
from netCDF4 import  Dataset # for using netCDF files
import numpy as np  # for scientific operators
from time import time  # for module1
from os import remove
import os as os
import pandas as pd  # conda install pandas
# for plotting
import seaborn as sns

from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap
import matplotlib.pyplot as plt
#import plotly.plotly as py

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from sherpa_auxiliaries import (create_emission_reduction_dict,
                                create_emission_dict, create_window,
                                read_progress_log, write_progress_log)
from sherpa_auxiliaries_epe import (gridint_toarray, write_nc, tiftogridgeneral,
                                    save_obj, load_obj, bm_precompute)
from sherpa_globals import (path_base_conc_cdf_test, path_salt_conc_cdf_test,
                            path_dust_conc_cdf_test, path_pop_mat_test,
                            path_mortbaseline, path_emission_cdf_test,
                            path_model_cdf_test, path_result_cdf_test)
from module7_custom import read_nuts_area

import module1 as shrp
#import module1old as shrp # at the moment using old model and module!


def blamematrix(perreduction, sources, path_results, precursors):
    """
    Calculates the blame matrix.
    Input:
        - perreduction: percentage reduction of the emissions
        - sources: list of countries
        - path_results
        - precursors
    Output:
        - df_avp: blamematrix dataframe concentration population weighted
                    normalized on delta emissions
        - df_avc: blamematrix dataframe concentration surface weighetd
                    normalized on delta emissions
        - df_c: blamematrix dataframe with only concentration changes

    Created on Thu Apr 20 16:26:42 2017
    @author: peduzem
    """

    # -------------------------------------------------------------------------
    # Calculating blamematrix
    # Warning: calculating emission delta for each precursor
    # only if it hasn't been already calc.
    # -------------------------------------------------------------------------
    if not os.path.exists('workdir\\dfemidelta.pkl'):
        dfemidelta = bm_precompute(perreduction, sources)
    else:
        print('WARNING: Loading previously calculated emissions deltas')
        dfemidelta = load_obj('dfemidelta')

    for prec in precursors:
        # Calculate dataframes with population, spatial average normalized by
        # emission reduction, and spatial average not normalized.
        df_avp, df_avc, df_c = bm_computation(dfemidelta, prec, sources)
        df_avp.to_csv(path_results + 'bm_avp_{}.csv'.format(prec))
        df_avc.to_csv(path_results + 'bm_avc_{}.csv'.format(prec))
        df_c.to_csv(path_results + 'bm_c_{}.csv'.format(prec))

def bm_computation(dfemidelta, prec, sources):
    """
    Calculates the blame matrix.
    Input:
        - dfemidelta: dataframe total reduction of precursor emission
        - sources: list of countries
        - precursor
    Output:
        -df: dataframes with the blame matrix
            df_avp: population weighted/emission reduction
            df_avc: averaged over the cells/emission reduction
            df_c: averaged over the cells

    Created on Thu Apr 20 16:26:42 2017
    @author: peduzem
    """
    targets = sources
    path_reduction_txt = 'workdir\\redbm_{}.txt'.format(prec)
    df_avp = pd.DataFrame(index=sources, columns=targets)
    df_avc = pd.DataFrame(index=sources, columns=targets)
    df_c = pd.DataFrame(index=sources, columns=targets)
    path_tiff = 'input/pop/7km_Qi_2010.tif'
    popall = tiftogridgeneral(path_tiff)
    path_surf_nc = 'input/JRC01.nc'
    # Area of each cell in the domain
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    # area_units = rootgrp.variables['surface'].units
    rootgrp.close()
    # -------------------------------------------------------------------------
    # Run module1 with the
    # emission reduction
    # -------------------------------------------------------------------------
    for source in sources:
        nc_redarea = 'workdir\\area_{}.nc'.format(source)
        # run module 1 with progress log
        # progresslog = 'input\\progress.log'
        output = 'blamematrix\\output_{}_{}\\'.format(prec, source)
        if not os.path.exists(output):
            os.makedirs(output)
            proglog_filename = path_result_cdf_test + 'proglog'
            write_progress_log(proglog_filename, 25, 2)
            start = time()
            shrp.module1(path_emission_cdf_test, nc_redarea,
                         path_reduction_txt, path_base_conc_cdf_test,
                         path_model_cdf_test, output)
#            shrp.module1(path_emission_cdf_test, nc_redarea,
#                         path_reduction_txt, path_model_cdf_test, output)
            stop = time()
            print('Module 1 for {}, {}, run time: {} sec.'.format(source, prec,
                                                              (stop-start)))
            remove(proglog_filename)
        else:
            print('WARNING: Loading previously calculated module1 results')
        # -------------------------------------------------------------------------

    for source in sources:
        print(prec)
        print(source)
        deltaconc = 'blamematrix\\output_{}_{}\\delta_concentration.nc'.format(prec, source)
        rootgrp = Dataset(deltaconc, mode='r')
        bc_pm25_conc = rootgrp.variables['delta_concentration'][:]
        # bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()
        for tar in targets:
            nc_area = 'workdir\\area_{}.nc'.format(tar)
            rootgrp = Dataset(nc_area, mode='r')
            area = rootgrp.variables['AREA'][:]
            # bc_pm25_units = rootgrp.variables['conc'].units
            rootgrp.close()
            conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))
            poptot = np.sum(popall * area / 100)
            df_avp.ix[tar][source] = (np.sum((1000 *popall * conc *
                                  area / 100) / poptot) / (
                                  dfemidelta.ix[source][prec] / 1000
                                  ))
            # Average over the surface (not over the cells)
            df_avc.ix[tar][source] = (np.sum((1000 * conc *
                                  area * surf/ 100) / np.sum(area * surf/ 100)) / (
                                  dfemidelta.ix[source][prec] / 1000
                                  ))
            df_c.ix[tar][source] = (np.sum((1000 *conc *
                                  area * surf / 100) / np.sum(area * surf/ 100)))
    return df_avp, df_avc, df_c

def bm_heatmap(df, name, precursor):
    if name =='avp':
        unit = '(ng/m3)/Gg'
        title = 'Normalized concentration change (population average)'
    elif name == 'avc':
        unit = '(ng/m3)/Gg'
        title = 'Normalized concentration change (concentration average)'
    elif name == 'c':
        unit = 'ng/m3'
        title = 'Concentration change'

    fig, ax = plt.subplots(figsize=(10, 10))
    cbar_ax = fig.add_axes([.905, 0.125, .03, 0.755])
    vmax = np.nanmax(df.values)/2
    vmin = np.nanmin(df.values)
    ax.set_title(('SHERPA, ' + title))
    ax = sns.heatmap(df, center=((vmin + vmax) / 2), vmax=vmax , annot=True, fmt='.1f', linewidths=.5, annot_kws={'size': 6}, cmap = cm.YlGnBu,
                     cbar_kws={'label': ('PM25 av. conc. change per {} em. change [{}]'.format(precursor,unit))}, ax=ax,  cbar_ax= cbar_ax)
    ax.set_xlabel('Source countries')  #
    ax.set_ylabel('Target countries')  #
    plt.savefig(path_results + 'SHERPAPM25{}_{}.png'.format(name, precursor),
                    dpi=300)
    sourcesemep=['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'GR', 'ES',
           'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
           'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'GB']
    df_diag = pd.Series(np.diag(df), index=[df.index, df.columns])
    if name == 'avc':
        if precursor == 'PPM':
            df_avc_fasst = pd.read_csv(path_input + 'EMEP.L00.PM25dPPM.csv',
                     index_col=[0])
        if precursor == 'NOx':
            df_avc_fasst = pd.read_csv(path_input + 'EMEP.L00.PM25dNOx.csv',
                     index_col=[0])
        if precursor == 'SOx':
            df_avc_fasst = pd.read_csv(path_input + 'EMEP.L00.PM25dSOx.csv',
                     index_col=[0])
        if precursor == 'NMVOC':
            df_avc_fasst = pd.read_csv(path_input + 'EMEP.L00.PM25dNMVOC.csv',
                     index_col=[0])
        if precursor == 'NH3':
            df_avc_fasst = pd.read_csv(path_input + 'EMEP.L00.PM25dNH3.csv',
                     index_col=[0])
        # Heat map
        df_avc_fasst_red=df_avc_fasst.ix[sourcesemep][sourcesemep]/0.15
        fig, ax = plt.subplots(figsize=(10, 10))
        cbar_ax = fig.add_axes([.905, 0.125, .03, 0.755])
        ax = sns.heatmap(df_avc_fasst_red, center=((vmin + vmax) / 2), vmax=vmax , annot=True, fmt='.1f', linewidths=.5, annot_kws={'size': 6}, cmap = cm.YlOrBr,
                 cbar_kws={'label': ('PM25 av. conc. change per {} em. change [{}]'.format(precursor,unit))}, ax=ax,  cbar_ax= cbar_ax)
        titlefasst = 'EMEP norm., ' + title
        ax.set_title(titlefasst)
        ax.set_xlabel('Source countries')  #
        ax.set_ylabel('Target countries')  #
        plt.savefig(path_results + 'FASST_{}.png'.format(precursor),dpi=300)
        # barplot on the diagonal
        df_fasst_diag = pd.Series(np.diag(df_avc_fasst_red), index=[df_avc_fasst_red.index, df_avc_fasst_red.columns])
        fig, ax = plt.subplots(figsize=(10,10))
        x_label = list(df.index.values)
        x=np.arange(1, (len(df.index)+1))
        y = df_diag.values
        z = df_fasst_diag.values
        ax = plt.subplot(111)
        ax.set_title('Comparison of concentration changes (with EMEP norm.)')
        rects1 = ax.bar(x, y,width=0.4,color='g',align='center')
        rects2 = ax.bar(x-0.4, z,width=0.4,color='r',align='center')
        ax.set_xlabel('Source/Target countries')
        ax.set_ylabel('PM25 av. conc. change per {} em. change [{}]'.format(precursor,unit))
        ax.set_xticklabels(x_label, rotation='vertical')
        ax.legend((rects1[0], rects2[0]), ('SHERPA-CHIMERE', 'EMEP norm.'))
        plt.xticks(np.arange(min(x), max(x)+1, 1.0))
        plt.savefig(path_results + 'diag{}_{}.png'.format(name,precursor), dpi=300)
        plt.show()

    if name == 'c':
        if precursor == 'PPM':
            df_c_emep = pd.read_excel(
                    path_results + 'input\\2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5', index_col=[0])
        if precursor == 'NOx':
             df_c_emep = pd.read_excel(
                    path_results + 'input\\2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5_NOx', index_col=[0])
        if precursor == 'SOx':
            df_c_emep = pd.read_excel(
                    path_results + 'input\\2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5_SOx', index_col=[0])
        if precursor == 'NMVOC':
            df_c_emep = pd.read_excel(
                    path_results + 'input\\2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5_NMVOC', index_col=[0])
        if precursor == 'NH3':
            df_c_emep = pd.read_excel(
                    path_results + 'input\\2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5_NH3', index_col=[0])
        # Heat map
        df_c_emep_red=df_c_emep.ix[sourcesemep][sourcesemep]
        fig, ax = plt.subplots(figsize=(10, 10))
        cbar_ax = fig.add_axes([.905, 0.125, .03, 0.755])
        ax = sns.heatmap(df_c_emep_red, center=((vmin + vmax) / 2), vmax=vmax , annot=True, fmt='.1f', linewidths=.5, annot_kws={'size': 6}, cmap = cm.YlOrBr,
                 cbar_kws={'label': ('PM25 av. conc. change per {} em. change [{}]'.format(precursor,unit))}, ax=ax,  cbar_ax= cbar_ax)
        titleemep = 'EMEP, ' + title
        ax.set_title(titleemep)
        ax.set_xlabel('Source countries')  #
        ax.set_ylabel('Target countries')  #
        plt.savefig(path_results + 'EMEP_{}.png'.format(precursor),dpi=300)
        # barplot on the diagonal
        df_emep_diag = pd.Series(np.diag(df_c_emep_red), index=[df_c_emep_red.index, df_c_emep_red.columns])
        fig, ax = plt.subplots(figsize=(10,10))
        x_label = list(df.index.values)
        x=np.arange(1, (len(df.index)+1))
        y = df_diag.values
        z = df_emep_diag.values
        ax = plt.subplot(111)
        ax.set_title('Comparison of concentration changes (with EMEP)')
        rects1 = ax.bar(x, y,width=0.4,color='g',align='center')
        rects2 = ax.bar(x-0.4, z,width=0.4,color='r',align='center')
        ax.set_xlabel('Source/Target countries')
        ax.set_ylabel('PM25 av. conc. change per {} em. change [{}]'.format(precursor,unit))
        ax.set_xticklabels(x_label, rotation='vertical')
        ax.legend((rects1[0], rects2[0]), ('SHERPA-CHIMERE', 'EMEP'))
        plt.xticks(np.arange(min(x), max(x)+1, 1.0))
        plt.savefig(path_results + 'diag{}_{}.png'.format(name,precursor), dpi=300)
        plt.show()

# -------------------------------------------------------------------------
# main program starts here
# -------------------------------------------------------------------------

if __name__ == '__main__':
    sources = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES',
               'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
               'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
#    for source in sources:
#        areacountry = gridint_toarray('NUTS_Lv0', 'parea', source)
#        path_areacountry_nc = 'workdir\\area_{}.nc'.format(source)
#        write_nc(areacountry, path_areacountry_nc, 'AREA', '%')

    perreduction = 15
    precursors = ['PPM', 'NOx', 'NMVOC', 'SOx','NH3']

    path_results = 'blamematrix\\'
    blamematrix(perreduction, sources, path_results, precursors)
## Normalized data:
    path_input = 'blamematrix\\input\\'
    sourcesemep=['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'GR', 'ES',
           'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
           'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'GB']
    for precu in precursors:

        if precu == 'PPM':
            df_c_emep = pd.read_excel(
                    path_input + '2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5', index_col=[0])
            df_c_emep_eu=df_c_emep.ix[sourcesemep][sourcesemep]
            df_e_emep = pd.read_excel(
                    path_input + 'emep_emissions.xlsx',
                    sheetname = 'PM2.5', index_col=[1])
            df_e_emep_eu = df_e_emep.ix[sourcesemep][2010]
        else:
             df_c_emep = pd.read_excel(
                    path_input + '2010_SRmatrices_R1Status2012AppC.xls',
                    sheetname = 'PM2.5_{}'.format(precu), index_col=[0])
             df_c_emep_eu=df_c_emep.ix[sourcesemep][sourcesemep]
             df_e_emep = pd.read_excel(
                    path_input + 'emep_emissions.xlsx',
                    sheetname = '{}'.format(precu), index_col=[1])
             df_e_emep_eu = df_e_emep.ix[sourcesemep][2010]
        df_emepl_eu = df_c_emep_eu / df_e_emep_eu
        df_emepl_eu.to_csv(path_input + 'EMEP.L00.PM25d{}.csv'.format(precu))

    names = ['avp', 'avc', 'c']
    for precu in precursors:
        print(precu)
        for name in names:
            print(name)
            df=pd.read_csv(path_results + 'bm_{}_{}.csv'.format(name, precu),
                     index_col=[0])
            bm_heatmap(df, name, precu)




        # Heat map
#        df_avc_emep_red=df_avc_emep.ix[sourcesemep][sourcesemep]
                # -------------------------------------------------------------------------

                        # uncomment to read results for debugging
##        df_avp = pd.read_csv(path_results + 'bm_avp_{}.csv'.format(precursor),
##                         index_col=[0])
##        df_avc = pd.read_csv(path_results + 'bm_avc_{}.csv'.format(precursor),
##                         index_col=[0])
##        df_c = pd.read_csv(path_results + 'bm_c_{}.csv'.format(precursor),
##                         index_col=[0])
#
#        # Preparing the data to make a plot, order on the diagaonal
#        df_avp_diag = pd.Series(np.diag(df_avp),
#                                index=[df_avp.index,
#                                       df_avp.columns]).sort_values()
#        ind_diag = list(df_avp_diag.index.values)
##        new_ind = [ind_diag[i][0] for i in np.arange(0, len(ind_diag))]
#        # TODO better figure
##        plotbm(df_avp.ix[new_ind][new_ind], 'avp', prec)
#
#        df_avc_diag = pd.Series(np.diag(df_avc),
#                                index=[df_avc.index,
#                                       df_avc.columns]).sort_values()
#        ind_diag = list(df_avc_diag.index.values)
##        new_ind = [ind_diag[i][0] for i in np.arange(0, len(ind_diag))]
##        plotbm(df_avc.ix[new_ind][new_ind], 'avc', prec)
#
#        df_c_diag = pd.Series(np.diag(df_c),
#                              index=[df_c.index,
#                              df_c.columns]).sort_values()
#        ind_diag = list(df_c_diag.index.values)
#        new_ind = [ind_diag[i][0] for i in np.arange(0, len(ind_diag))]
#        plotbm(df_c.ix[new_ind][new_ind], 'c', prec)

#    def plotbm(df, name, prec):
#        """
#        Plots the data frame for each precurosr as a 3D surface
#        """
#
##        fig, ax = plt.subplots(figsize=(10, 10))
##    #    sns.palplot(sns.color_palette("cubehelix", 80))
##
###        vmax = np.nanmax(df.values)/2
##        ax = sns.heatmap( df_avc, vmax=200, annot=True, fmt='.1f', linewidths=.5, annot_kws={'size': 6},
##                         cbar_kws={'label': 'PM25 av. conc. change per PPM em. change [(ng/m3)/Gg]'}, ax=ax)
##    #    cmap = 'plasma'
##
##        ax.set_xlabel('Source countries')  #
##        ax.set_ylabel('Target countries')  #
##        plt.savefig(path_results + 'SHERPAPM25_{}.png'.format(precursor),
##                        dpi=300)
##        plt.show()
#
#        fig = plt.figure(figsize=(10, 8))
#        ax = fig.gca(projection='3d')
#        fontsize = 6
#        # Make data.
#        ylabel = list(df.index.values) # TARGETS
#        xlabel = list(df.columns.values) # SOURCES
#        x = np.arange(len(xlabel))
#        y = np.arange(len(ylabel))
#        x, y = np.meshgrid(x, y)
#        # Plot the surface.
#        surf = ax.plot_surface(x, y, df.values, cmap=cm.coolwarm, antialiased=True)
##                               linewidth=0, vmax=200) #,
#        # Set ticks labels for x-axis
#        plt.xticks(np.arange(min(x[0]), max(x[0])+1, 1.0))
#        plt.yticks(np.arange(min(x[0]), max(x[0])+1, 1.0))
#        ax.set_xlabel('Source countries')  # TODO check
#        ax.set_ylabel('Target countries')  # TODO check
#        if name == 'avc':
#            ax.set_zlabel(r'Average conc change [(ng/m3)/Gg]')
#        if name == 'avp':
#            ax.set_zlabel(r'Population weighted change [((ng/m3)/Gg]')
#        if name == 'c':
#            ax.set_zlabel(r'Not-normalized change [(ng/m3)]')
#        ax.set_xticklabels(xlabel, rotation='vertical')
#        ax.set_yticklabels(ylabel, rotation='vertical')
#        # set fontsize
#        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label,
#                     ax.zaxis.label] + ax.get_xticklabels() +
#                     ax.get_yticklabels() + ax.get_zticklabels()):
#            item.set_fontsize(fontsize)
#        # Add a color bar which maps values to colors.
#        cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
#        if name == 'c':
#            cbar.set_label(r'Average conc change [(ng/m3)]',
#                           fontsize=fontsize)
#        else:
#            cbar.set_label(r'Average conc change [(ng/m3)/Gg]',
#                           fontsize=fontsize)
#        cbar.ax.tick_params(labelsize=fontsize)
#        plt.savefig(path_results + name + 'bm_{}.png'.format(prec),
#                    dpi=300)
#        plt.show()

#        ax = sns.heatmap(df)
