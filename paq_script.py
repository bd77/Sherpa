# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 09:39:24 2017

Script for the Partnership of Air Quality
It provides the impact (yll or mortality) on the Functional Urban Area of
different cities given a reduction of emission at different levels: the city,
the FUA, the country and the whole of Europe.

@author: peduzem
"""

from os import remove  # for module1
from time import time  # for module1

import matplotlib.pyplot as plt
from netCDF4 import Dataset  # for using netCDF files
import numpy as np
import pandas as pd


from healthia import calc_impacts
import module1old as shrp # at the moment using old model and module!
from sherpaeco import save_obj, load_obj
from sherpa_auxiliaries import (create_emission_reduction_dict,
                                create_emission_dict, create_window,
                                read_progress_log, write_progress_log)
from sherpa_auxiliaries_epe import (write_nc, gridint_toarray)
from sherpa_globals import (path_area_cdf_test, path_model_cdf_test,
                            path_emission_cdf_test, path_result_cdf_test,
                            path_base_conc_cdf_test)

if __name__ == '__main__':

    # -------------------------------------------------------------------------
    # Preparing the input data
    # -------------------------------------------------------------------------
    print('preparing input data...')
    # Open model file, I just need it to load lons and lats etc.
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(lon_array)
    n_lat = len(lat_array)
    rootgrp.close()

    # Area of each cell in the domain
    nc_area = 'input\\JRC01.nc'
    rootgrp = Dataset(nc_area, mode='r')
    area = rootgrp.variables['surface'][:]
    area_units = rootgrp.variables['surface'].units
    rootgrp.close()

#    # Dust and sald PM25 conentrations to be removed for HIA
#    nc_dust = 'input\\pDUST-pSALT\\pDUST-25-basecase.nc'
#    rootgrp = Dataset(nc_dust, mode='r')
#    pDUST25 = rootgrp.variables['pDUST-25'][:]
#    pDUST25_units = rootgrp.variables['pDUST-25'].units
#    rootgrp.close()
#
#    nc_salt = 'input\\pDUST-pSALT\\pSALT-25-basecase.nc'
#    rootgrp = Dataset(nc_salt, mode='r')
#    pSALT25 = rootgrp.variables['pSALT-25'][:]
#    pSALT25_units = rootgrp.variables['pSALT-25'].units
#    rootgrp.close()

#    rootgrp = Dataset(path_base_conc_cdf_test, mode='r')
#    bc_conc = rootgrp.variables['conc'][:]
#    bc_conc_units = rootgrp.variables['conc'].units
#    rootgrp.close()
#
#    nc_antrconc='workdir\\antrconc.nc'
#    fh = Dataset(nc_antrconc, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude', len(lat_array))
#    fh.createDimension('longitude', len(lon_array))
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    antrconc = fh.createVariable('conc', 'f8', ('latitude', 'longitude'))
#    fh.variables['conc'].units =  bc_conc_units
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    antrconc[:] = bc_conc - pSALT25 - pDUST25
#    fh.close()

    # Selection of the EU28 domain
    eu28_codes = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES',
                  'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
                  'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']

    area_eu28 = np.zeros((1, len(lat_array), len(lon_array)))
    for eucode in eu28_codes:
        area_eu28 = area_eu28 + (gridint_toarray(
                                 'NUTS_Lv0', 'parea', eucode))
#    area_eu28 = np.ones((1, len(lat_array), len(lon_array)))*100

    path_area_eu = 'workdir\\redarea_eu'
    write_nc(area_eu28, path_area_eu, 'AREA', '%')

    # specify the path to the reduction
    path_reduction_txt = 'input\\paqred.txt'

    # list of cities in the Partenership of air quality
    # Cities:
    cities = ['Helsinki', 'Milano', 'Constanta', 'Grad Zagreb',
              'Duisburg', 'London', 'Ostrava', 'Utrecht']
    # Corresponding codes for FUA (target areas)
    # NB: for Duisburg the FUA name is Ruhrgebiet
    tarcodes = ['FI001L2', 'IT002L2', 'RO501L1', 'HR001L2', 'DE038L1',
                'UK001L2', 'CZ003L1', 'NL004L2']
    # Corresponding codes for great cities (reduction areas)
    redgcodes = ['FI001K2', 'IT002K1', 'RO501C1', 'HR001C1', 'DE501C1',
                 'UK001K2', 'CZ003C1', 'NL004C1']
    redfuacodes = tarcodes  # Corresponding codes for FUA
    redcocodes = ['FI', 'IT', 'RO', 'HR', 'DE', 'UK', 'CZ', 'NL']
    redlevels = ['GCITY_CODE', 'FUA_CODE', 'NUTS_Lv0']

        # list of cities in the Partenership of air quality
#    # Cities:
#    cities = ['Helsinki'] #, 'Milano', 'Constanta', 'Grad Zagreb',
#              #'Duisburg', 'London', 'Ostrava', 'Utrecht']
#    # Corresponding codes for FUA (target areas)
#    # NB: for Duisburg the FUA name is Ruhrgebiet
#    tarcodes = ['FI001L2'] #, 'IT002L2', 'RO501L1', 'HR001L2', 'DE038L1',
#                #'UK001L2', 'CZ003L1', 'NL004L2']
#    # Corresponding codes for great cities (reduction areas)
#    redgcodes = ['FI001K2'] #,, 'IT002K1', 'RO501C1', 'HR001C1', 'DE501C1',
#                 #'UK001K2'] #,, 'CZ003C1', 'NL004C1']
#    redfuacodes = tarcodes  # Corresponding codes for FUA
#    redcocodes = ['FI'] #,, 'IT', 'RO', 'HR', 'DE', 'UK', 'CZ', 'NL']
#    redlevels = ['GCITY_CODE', 'FUA_CODE', 'NUTS_Lv0']

    print('starting')

    # DataFrame with the city names (cities) as columns and as
    # indeces the reduction levels.
    df = pd.DataFrame(([redgcodes, redfuacodes, redcocodes]),
                      index=redlevels, columns=cities)

    # Creating a list of tuples to access the FUA target areas
    len_tar = len(df.loc['FUA_CODE'])
    pareas = ['parea']*len_tar
    tartuple = list(zip(['FUA_CODE']*len_tar, pareas, tarcodes))

    nc_tararea = {}  # path to store the target area netcdf file
    for cityind in np.arange(len(cities)):
        # path to store the target areas
        nc_tararea[cityind] = (
         'workdir\\tararea_{}.nc'.format(tartuple[cityind][2]))
        # prepare the target areas, from the grid interesect using the
        # tartuple (target tuples)
        area_sel = gridint_toarray(tartuple[cityind][0], tartuple[cityind][1],
                                   tartuple[cityind][2])
        write_nc(area_sel, nc_tararea[cityind], 'AREA', '%')

    # initialize dictionaries to store:
    nc_redarea = {}  # path of netcdf where there are the reductions
    deltayll_reg = {}  # the delta yll achieved in the target area
    deltayll_tot = {}  # the delta yll achieved in total
    delta_mort_tot = {}  # the delta deaths achieved in total
    delta_mort_reg = {}  # the delta deaths achieved in the target area
    deltayll_spec_reg = {}  # the array delta yll achieved per cell

    # Preapre all the areas of reduction
    redtuple = {}
    for i in redlevels:
        # Creating a list of tuples to access the reduction areas
        redtuple[i] = list(zip([i]*len_tar, pareas, df.loc[i]))
        # path to store the selected area
        for cityind in np.arange(len(cities)):
            nc_redarea = 'workdir\\redarea{}_{}.nc'.format(
                    i, redtuple[i][cityind][2])
            redarea = gridint_toarray(redtuple[i][cityind][0], redtuple[i][cityind][1],
                            redtuple[i][cityind][2])
            write_nc(redarea, nc_redarea, 'AREA', '%')

    # -------------------------------------------------------------------------
    # Running module1 and calculating impacts
    # -------------------------------------------------------------------------
    for i in redlevels:
        deltayll_reg[i] = {}
        deltayll_tot[i] = {}
        delta_mort_tot[i] = {}
        delta_mort_reg[i] = {}
        deltayll_spec_reg[i] = {}
        for cityind in np.arange(len(cities)):
            # read reduction areas netcdf fiels
            nc_redarea = 'workdir\\redarea{}_{}.nc'.format(
                    i, redtuple[i][cityind][2])
            # run module 1 with progress log
            progresslog = 'input\\progress.log'
            output = 'output\\'
            proglog_filename = path_result_cdf_test + 'proglog'
            write_progress_log(proglog_filename, 25, 2)
            start = time()
#            shrp.module1(path_emission_cdf_test, nc_redarea,
#                         path_reduction_txt, path_base_conc_cdf_test, path_model_cdf_test, output)
            shrp.module1(path_emission_cdf_test, nc_redarea,
                         path_reduction_txt, path_model_cdf_test, output)
            stop = time()
            print('Module 1 run time: %s sec.' % (stop-start))
            remove(proglog_filename)
            deltaconc = 'output\\delta_concentration.nc'

            # Evaluation of the impacts
            (deltayll_reg[i][tartuple[cityind][2]],
             delta_mort_reg[i][tartuple[cityind][2]],
             deltayll_spec_reg[i][tartuple[cityind][2]]) = calc_impacts(
                     deltaconc, nc_tararea[cityind], tartuple[cityind][2][0:2])

    # impact of basecase (should maybe eslcude natural?)
    deltayll_reg_bc = {}
    deltayll_tot_bc = {}
    delta_mort_tot_bc = {}
    delta_mort_reg_bc = {}
    deltayll_spec_reg_bc = {}

    # impact of reductions at european level
    deltayll_reg_eu = {}
    deltayll_tot_eu = {}
    delta_mort_tot_eu = {}
    delta_mort_reg_eu = {}
    deltayll_spec_reg_eu = {}

    for cityind in np.arange(len(cities)):
        (deltayll_reg_bc[cityind], delta_mort_reg_bc[cityind],
         deltayll_spec_reg_bc[cityind]) = calc_impacts(
                 path_base_conc_cdf_test, nc_tararea[cityind], tartuple[cityind][2][0:2], spec='a')
        # run module 1 with progress log
        progresslog = 'input\\progress.log'
        output = 'output\\'
        proglog_filename = path_result_cdf_test + 'proglog'
        write_progress_log(proglog_filename, 25, 2)
        start = time()
#        shrp.module1(path_emission_cdf_test, path_area_eu,
#                     path_reduction_txt, path_base_conc_cdf_test, path_model_cdf_test, output)
        shrp.module1(path_emission_cdf_test, path_area_eu,
                     path_reduction_txt, path_model_cdf_test, output)
        stop = time()
        print('Module 1 run time: %s sec.' % (stop-start))
        remove(proglog_filename)
        deltaconc = 'output\\delta_concentration.nc'
        # Evaluation of the impacts
        (deltayll_reg_eu[cityind], delta_mort_reg_eu[cityind],
         deltayll_spec_reg_eu[cityind]) = calc_impacts(
                 deltaconc, nc_tararea[cityind], tartuple[cityind][2][0:2], spec='d')

    # -------------------------------------------------------------------------
    # Display results
    # -------------------------------------------------------------------------
    print('Displaying results in YLL')
    for cityind in np.arange(len(cities)):
#        cityind = 0
        yllred = [deltayll_reg[i][tartuple[cityind][2]] for i in redlevels]
        yllreddelta = (deltayll_reg_bc[cityind][1] - [yllred[i][1] for i in np.arange(0,3)])

        yllredplot = [
                deltayll_reg_bc[cityind][1],
                deltayll_reg_bc[cityind][1] - deltayll_reg_eu[cityind][1]]

        yllredplot[1:1] = yllreddelta
        N = len(yllredplot)
        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind[0], yllredplot[0], width, color='r')
        rects1 = ax.bar(ind[1:], yllredplot[1:], width, color='g')

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Years of life loss [YLL]')
        ax.set_title(
         'Years of life loss in the FUA of {}, year 2010'.format(cities[cityind]))
        ax.set_xticks(ind)
        ax.set_xticklabels((['BC', 'city', 'FUA', 'country', 'Europe']))
        ax.legend(
                ('Base Case', '30% reduction'),
                bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(
         'output\\yll{}.png'.format(
                cities[cityind]), bbox_inches='tight', dpi=300)

    print('Displaying results in YLL per 100000 ppl')
    for cityind in np.arange(len(cities)):
#        cityind = 0
        yllreds = [
                deltayll_spec_reg[i][tartuple[cityind][2]] for i in redlevels]
        yllreddeltas = (deltayll_spec_reg_bc[cityind][1]-[yllreds[i][1] for i in np.arange(0,3)])
#
        yllredplots = [
                deltayll_spec_reg_bc[cityind][1],
                deltayll_spec_reg_bc[cityind][1] - deltayll_spec_reg_eu[cityind][1]]
#
        yllredplots[1:1] = yllreddeltas
        N = len(yllredplots)
        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars
#
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind[0], yllredplots[0], width, color='r')
        rects1 = ax.bar(ind[1:], yllredplots[1:], width, color='g')

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Specifi years of life loss [YLL per 100k ppl]')
        ax.set_title(
         'Years of life loss per 100k ppl (FUA, year 2010, 30% reductions)'.format(
                 cities[cityind]))
        ax.set_xticks(ind)
        ax.set_xticklabels((['BC', 'city', 'FUA', 'country', 'Europe']))
        ax.legend(
                ('Base Case', '30% reduction'),
                bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(
         'output\\yll_spec{}.png'.format(
                 cities[cityind]), bbox_inches='tight', dpi=300)
        plt.show()

    print('Displaying results in YLL per 100000 ppl, in one figure')
    for cityind in np.arange(len(cities)):
#        cityind = 0
        yllreds = [
                deltayll_spec_reg[i][tartuple[cityind][2]] for i in redlevels]
        yllreddeltas = (deltayll_spec_reg_bc[cityind][1]-[yllreds[i][1] for i in np.arange(0,3)])
#
        yllredplots = [
                deltayll_spec_reg_bc[cityind][1],
                deltayll_spec_reg_bc[cityind][1] - deltayll_spec_reg_eu[cityind][1]]
#
        yllredplots[1:1] = yllreddeltas
        N = len(yllredplots)
        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        ax = plt.subplot(111)
        ax.set_ylabel('Specific years of life loss [YLL per 100k ppl]')
        ax.set_title('Years of life loss per 100k ppl (FUA, year 2010, 30% reductions)')
        ax.set_xticks(ind)
        ax.set_xticklabels((['BC', 'city', 'FUA', 'country', 'Europe']))
        plt.plot(ind, yllredplots, label='{}'.format(cities[cityind]), marker='o')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig(
         'output\\yll_spec_comp.png', bbox_inches='tight', dpi=300)

    print('Displaying results in YLL, relative')
    for cityind in np.arange(len(cities)):
        yllred = [deltayll_reg[i][tartuple[cityind][2]] for i in redlevels]
        yllreddelta = (deltayll_reg_bc[cityind][1] - [yllred[i][1] for i in np.arange(0,3)])
        yllreddeltarel= (yllreddelta /
                       deltayll_reg_bc[cityind][1])

        yllredplotrel = [
                deltayll_reg_bc[cityind][1]/deltayll_reg_bc[cityind][1],
                (deltayll_reg_bc[cityind][1] -
                 deltayll_reg_eu[cityind][1]) / deltayll_reg_bc[cityind][1]]

        yllredplotrel[1:1] = yllreddeltarel

        N = len(yllredplotrel)
        ind = np.arange(N)

#        fig, ax = plt.subplots()

        ax = plt.subplot(111)
#        ax = plt.Axes
        ax.set_ylabel('Relative years of life loss [YLL/YLL]')
        ax.set_title('Relative years of life loss (FUA, year 2010, 30% reductions)')
        ax.set_xticks(ind)
        ax.set_xticklabels((['BC', 'city', 'FUA', 'country', 'Europe']))
        plt.plot(ind, yllredplotrel, label='{}'.format(cities[cityind]), marker='o')
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig(
         'output\\yll_relative.png', bbox_inches='tight', dpi=300)


    print('Displaying results in Mortality')

    for cityind in np.arange(len(cities)):
        mortred = [
                delta_mort_reg[i][tartuple[cityind][2]] for i in redlevels]
        mortreddelta = (delta_mort_reg_bc[cityind][1]-[mortred[i][1] for i in np.arange(0,3)])

        mortredplot = [
                delta_mort_reg_bc[cityind][1],
                delta_mort_reg_bc[cityind][1] - delta_mort_reg_eu[cityind][1]]

        mortredplot[1:1] = mortreddelta
        # lower bound
#        lmortred = [
#                delta_mort_reg[i][tartuple[cityind][2]][0] for i in redlevels]
        lmortreddelta = (delta_mort_reg_bc[cityind][0]-[mortred[i][0] for i in np.arange(0,3)])

        lmortredplot = [
                delta_mort_reg_bc[cityind][0],
                delta_mort_reg_bc[cityind][0] - delta_mort_reg_eu[cityind][0]]

        lmortredplot[1:1] = lmortreddelta

        # upper bound
#        umortred = [
#                delta_mort_reg[i][tartuple[cityind][2]][2] for i in redlevels]
        umortreddelta = (delta_mort_reg_bc[cityind][2]-[mortred[i][2] for i in np.arange(0,3)])

        umortredplot = [
                delta_mort_reg_bc[cityind][2],
                delta_mort_reg_bc[cityind][2] - delta_mort_reg_eu[cityind][2]]

        umortredplot[1:1] = umortreddelta
        # confidence interval
        ci = np.asarray(umortredplot)-np.asarray(lmortredplot)

        N = len(mortredplot)
        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind[0], mortredplot[0], width, color='r') # yerr = ci[0],
        rects1 = ax.bar(ind[1:], mortredplot[1:], width, color='g') # yerr = ci[1:],

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Mortality [number of deaths]')
        ax.set_title(
         'Mortality in the FUA of {}'.format(cities[cityind]))
        ax.set_xticks(ind)
        ax.set_xticklabels((['BC', 'city', 'FUA', 'country', 'Europe']))
        ax.legend(
                ('Base Case', '30% reduction'),
                bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(
         'output\\mort{}.png'.format(cities[cityind]), bbox_inches='tight', dpi=300)

    ## quick checks
#    nc_area='workdir\\redareaNUTS_Lv0_IT.nc'
#    nc_redarea='workdir\\redareaNUTS_Lv0_IT.nc'
#    nc_redarea = 'workdir\\redareaFUA_CODE_UK001L2.nc'
#    nc_redarea = 'workdir\\redareaGCITY_CODE_UK001K2.nc'
#    concncdf = nc_antrconc
#    (deltayll_reggl, deltayll_totgl,
#     delta_mort_totgl, delta_mort_reggl,
#     delta_yllgl) = calc_impacts(
#             path_base_conc_cdf_test, area_eu, spec='a')