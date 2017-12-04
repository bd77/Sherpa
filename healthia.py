# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017

Collection of functions to evaluate the helath impact of AQ
(TODO: NO2, and PM10... at the moment only PM2.5 is considered)

@author: peduzem
"""
from netCDF4 import Dataset  # for using netCDF files
import numpy as np
# from osgeo import gdal, ogr, osr  # conda install -c conda-forge gdal
import pandas as pd
# import scipy.io as sio
import os as os

from sherpa_auxiliaries_epe import gridint_toarray, write_nc, tiftogridgeneral


def calc_impacts(path_conc_nc, path_areasel_nc, path_tiff, path_mortbaseline,
                 path_dust_conc_cdf_test, path_salt_conc_cdf_test,
                 code, spec='d',
                 approx='e', miller=True,
                 std_life_exp=70):
    """
    Health impacts of PM2.5 - WHO-Europe method (HRAPIE reccomendations)

    Concentration-Response-Function taken from the software AirQ+ (WHO)

    YLLs calculated considering the distribution of mortality and population
    by age and by country, where data is not available in the ICD-10 format
    YLLs can be calculated as a function of life expectancy for PM2.5 or
    considering average EU values.

    NB: Hereafter a positive delta means a reduction!

    Input:
        - path_conc_nc: path to the file  with the concentrations
           (output of module1) or absolute concentration
        - path_areasel_nc: path to netcdf file of area of interest
           (nc file with 100 in the cells of the area of interest)
        - path_tiff: .tif population file (i.e. LUISA)
        - path_mortbaseline: excel file with the data of the baseline
        - path_dust_conc_cdf_test: concentration of natural dust
        - path_salt_conc_cdf_test: concentration of sea salt
        - code: country code, FUA code or city code (as defined in the
           grid interesect files):
        - spec='d' for 'delta' concentration or 'a' for 'absolute'
        - std_life_exp=70: standard life expectancy to calculate the years of
           life loss, 70 is the value used by WHO in
           http://data.euro.who.int/dmdb/
        - approx = 'e' for exponential and 'l' for linear
        - miller = True to use the aproxmation of YLL as described in Holland
          (I could not find the original miller reference) for countries
          for which I don't have data.

    Output:
        - deltayll_reg: total years of life loss in the selected area
        - delta_mort_reg: total deaths in the selected area
        - deltayll_spec_reg: YLL per 100000 ppl in the selected area

    Required files:
        - path_mortbaseline: excel file with the data of the baseline
        - path_tiff: .tif file with the population distribution (LUISA)
        - grid intersect files
        Country codes for teh grid intesect are:
             'CY', 'ES', 'HU', 'ME', 'NL', 'SI', 'AT', 'BE', 'BG', 'CH', 'CZ',
             'DE', 'DK', 'EE', 'EL', 'FI', 'FR', 'HR', 'IE', 'IT',
             'LT', 'LU', 'LV', 'MK', 'MT', 'NO', 'PL', 'PT', 'RO', 'SE', 'SK',
             'TR', 'UK'.
        However some of these countries have no or incomplete population data:
                 'CH', 'MK', 'ME', 'NO', 'TR'. 'SE', 'NO',

    BIBLIOGRAPHY:

    - Methodology:
        World Health Organization Europe, 2013. Health risks of air pollution
        in Europe - HRAPIE project - Recommendations for concentration–response
        functions for cost–benefit analysis of particulate matter, ozone and
        nitrogen dioxide, Copenhagen Ø, Denmark.

        Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios
        for the EU Clean Air Package Version. Version 2

    - Concentration-response function
        World Health Organization Europe, 2017. AirQ+: software tool for health
        risk assessment of air pollution.

    - Data for baseline population (path_mortbaseline):
        ICD codes: ICD-10: A00-B99,C00-D48,D50-D89,E00-E88,F01-F99,G00-G98,
        H00-H59,H60-H93,I00-I99,J00-J98,K00-K92,L00-L98,M00-M99,N00-N98,
        O00-O99,P00-P96,Q00-Q99,R00-R99
        Age: '30 - 85 +'
        Sex: Both
        http://data.euro.who.int/dmdb/ [Accessed December 13, 2016].

        Life expectancy at bith for the year 2009 (most recent with the most
        values) from:
        http://data.euro.who.int/hfamdb/ [Accessed Aprli 12, 2017].

    @author: peduzem
    """

    # -----------------------------------------------------------------------------
    # BASELINE POPULATION DATA
    # Load population file from LUISA
    popall = tiftogridgeneral(path_tiff)
    # remove all negative and -inf values (not sure why they are there)
    popall = np.where(np.isinf(popall), 0, popall)
    popall = np.where(popall < 0, 0, popall)
#    write_nc(popall, 'input\\pop\\pop_nc.nc', 'ppl', '-')
    # Build Pandas dataframe with the baseline population data
    # distributed by age:
    # Baseline morality
    df_mort = pd.read_excel(path_mortbaseline,
                            sheetname='mort', skiprows=4,  header=[1],
                            index_col=[0])
    # Baseline population
    df_pop = pd.read_excel(path_mortbaseline,
                           sheetname='pop', skiprows=4,  header=[1],
                           index_col=[0])

    # rename rows : country names with corresponding country codes
    # rename columns: age intervals with mid point age
    df_cn = pd.read_excel(path_mortbaseline,
                          sheetname='map_country',
                          index_col=[0])
    df_y = pd.read_excel(path_mortbaseline,
                         sheetname='map_years',
                         index_col=[0])
    dic_cn = df_cn.to_dict()
    dic_y = df_y.to_dict()
    df_mort.rename(index=dic_cn['country_code'], inplace=True)
    df_pop.rename(index=dic_cn['country_code'], inplace=True)
    df_mort.rename(columns=dic_y['mid_year'], inplace=True)
    df_pop.rename(columns=dic_y['mid_year'], inplace=True)

    #  Baseline population and deaths older than 30 and younger
    #  than the std_life_exp
    col_list = list(df_mort)
    cols = [
            x for x in np.asarray(col_list[1:])
            if (x < std_life_exp and x > 30)]
    df_mort['sum_sel'] = df_mort[cols].sum(axis=1)
    df_pop['sum_sel'] = df_pop[cols].sum(axis=1)
    #  Life years lost correponding to each mid point age
    lyl = std_life_exp - np.asarray(cols)

    def lylcal(x):
        """
        Life years loss caclculation: average life year lost for each death
        older than 30 and younger than the std_life_exp, we can use this
        because the CRF does not depend on age
        """
        return np.sum(x[cols]*lyl/x['sum_sel'])

    df_lyl = df_mort.apply(lylcal, axis=1)

    # Death rate: all mortality and all population above 30
    # is accounted for
    cols30p = [x for x in np.asarray(col_list[1:]) if x > 30]
    drate = df_mort[cols30p].sum(axis=1)/df_pop[cols30p].sum(axis=1)

    # Ratio between the population over 30 and the total population
    # HP: The population over 30 has the same distribution as all population
    # (as provided by LUISA)
    p30_ptot = (df_pop[np.asarray(cols30p)].sum(axis=1) /
                df_pop[np.asarray(col_list[1:])].sum(axis=1))
#    ptot = df_pop[np.asarray(col_list[1:])].sum(axis=1)
    # Calculting average values to use if the country selected is not in
    # the database. Average values calculating considering EU28 countries
    # which are also in the WHO database

    eu34_codes =['AT', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL',  
                 'ES', 'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LI', 'LT', 'LU', 
                 'LV', 'ME', 'MK', 'MT', 'NL', 'NO', 'PL', 'PT', 'RO', 'SE',
                 'SI', 'SK', 'TR', 'UK']
    
    # Find coutnries that are both in the EU and in the mortality database
    intersect = list(
            set.intersection(set(eu34_codes), set(list(df_mort.index))))
    # Calculate averages for the countries in both lists:
    m_drate = drate.ix[intersect].mean(axis=0)
    m_p30_ptot = p30_ptot.ix[intersect].mean(axis=0)
    m_lyl = df_lyl.ix[intersect].mean(axis=0)
#    ptot_eu = df_pop.ix[intersect][col_list[1:]].sum(axis=1)
    # total population by country in both list

# -----------------------------------------------------------------------------
# PM2.5 CONCENTRATION

    if spec == 'a':
        # base case PM25 conentration
        rootgrp = Dataset(path_conc_nc, mode='r')
        bc_pm25_conc = rootgrp.variables['conc'][:]
#        bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()

        # Dust PM25 conentration to be removed for HIA
        rootgrp = Dataset(path_dust_conc_cdf_test, mode='r')
        pDUST25 = rootgrp.variables['pDUST-25'][:]
#        pDUST25_units = rootgrp.variables['pDUST-25'].units
        rootgrp.close()

        # Salt PM25 conentration to be removed for HIA
        rootgrp = Dataset(path_salt_conc_cdf_test, mode='r')
        pSALT25 = rootgrp.variables['pSALT-25'][:]
#        pSALT25_units = rootgrp.variables['pSALT-25'].units
        rootgrp.close()

        # Antropogenic concentration
        pm25_conc = bc_pm25_conc - pSALT25 - pDUST25

    elif spec == 'd':
        fh_deltapm25 = Dataset(path_conc_nc, mode='r')
        pm25_conc = fh_deltapm25.variables['delta_concentration'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
        fh_deltapm25.close()

# -----------------------------------------------------------------------------
    #  Get information sepcific to the country of the are of interest
    country = code[0:2]

    # area_co = gridint_toarray(
    #       'NUTS_Lv0', 'parea', country)
    # path_areco_nc = 'workdir/areaco{}.nc'.format(country)
    # write_nc(area_co, path_areco_nc, 'AREA', '%')
    # total population in the country according to LUISA
    pop30plus = np.zeros(np.shape(popall))
#    popallco = np.sum(popall * area_co/100)
    if country in list(df_mort.index):
        # Scaling factor to respect country totals as provided by WHO
        # in each grid cell of the area of interest
        # scalefact = ptot_eu.ix[country] / popallco
        # Scaling to only the population over 30 in each grid cell
        # warning: this is valid only in the country
        pop30plus = pop30plus + np.array((popall *
                                          p30_ptot.ix[country]))
    else:
        # Scaling to only the population over 30 in each grid cell
        # warning: this is valid only in the country
        pop30plus = pop30plus + np.array(popall * m_p30_ptot)

    # CONCENTRATION RESPONSE FUNCTION:
    # considering bounds for 95% CI
    if approx == 'l':
        # linear approximation
#        mrr = 1.06  # 'middle' value
#        lrr = 1.04  # lower bound
#        hrr = 1.083  # higher bound
#        # crf = 0.006
#        af = np.asarray([(lrr-1)/lrr, (mrr-1)/mrr, (hrr-1)/hrr]) / 10
#        pt = len(af)
#        if country in list(df_mort.index):
#            delta_mort = np.where(np.isnan(pm25_conc), 0, (
#                    [af[i]*pm25_conc*pop30plus*drate.ix[country]
#                     for i in range(len(af))]))
#        else:
#            delta_mort = np.where(np.isnan(pm25_conc), 0, (
#                    [af[i]*pm25_conc*pop30plus*m_drate
#                     for i in range(len(af))]))
        mrr = 0.06  # 'middle' value
        lrr = 0.04  # lower bound
        hrr = 0.083  # higher bound
        # crf = 0.006
        af = [lrr, mrr, hrr]
        pt = len(af)
        if country in list(df_mort.index):
            delta_mort = np.where(np.isnan(pm25_conc), 0, (
                    [(af[i]*pm25_conc/10)/(1+af[i]*pm25_conc/10)*pop30plus*drate.ix[country]
                     for i in range(len(af))]))
        else:
            delta_mort = np.where(np.isnan(pm25_conc), 0, (
                    [(af[i]*pm25_conc/10)/(1+af[i]*pm25_conc/10)*pop30plus*m_drate
                     for i in range(len(af))]))

    elif approx == 'e':
        # Taken from AirQ+ (WHO)
        # cutoff = 10 # microg/m3 # Taken from AirQ+ (WHO)
        # but not mentioned in the guidelines so not considered here.
        mbeta = 0.006015392281974714  # 'middle' value
        lbeta = 0.003922071315328133  # lower bound
        hbeta = 0.007973496801885352  # higher bound
        beta = [lbeta, mbeta, hbeta]
        pt = len(beta)
        if country in list(df_mort.index):
            delta_mort = np.where(
                    np.isnan(pm25_conc), 0, (
                            [(1-(np.exp(-beta[i]*pm25_conc))) *
                             pop30plus * drate.ix[country]
                             for i in range(len(beta))]))
        else:
            delta_mort = np.where(
                    np.isnan(pm25_conc), 0, (
                            [(1-(np.exp(-beta[i]*pm25_conc))) *
                             pop30plus * m_drate
                             for i in range(len(beta))]))

# -----------------------------------------------------------------------------
    # ESTIMATE OF THE YLL
    if country in list(df_mort.index):
        # years of life lost
        delta_yll = np.where(np.isnan(pm25_conc), 0,
                             [delta_mort[i] * df_lyl.ix[country]
                             for i in range(pt)])

    else:
        # if there is no data use the average values
        delta_yll = np.where(np.isnan(pm25_conc), 0,
                             [delta_mort[i] * m_lyl
                              for i in range(pt)])
        # if miller is true use the function in Holland by Miller
        if miller is True:
            df_le = pd.read_excel(path_mortbaseline,
                                  sheetname='le', skiprows=0, header=[0],
                                  parse_cols=[6, 8])
            # life years gained per 100000 ppl for a unit concentration change,
            # We don't have the confidence interval so the same value is
            # repeated 3 times
            df_le.index = df_le['co_code']
            lyg = np.exp(8.161-(0.04478*df_le[2009].ix[country]))
            delta_yll = np.where(np.isnan(pm25_conc), 0, [pm25_conc *
                                 lyg/100000*pop30plus for i in range(pt)])

    # PRODUCE RESULTS for the area of interest
    rootgrp = Dataset(path_areasel_nc, mode='r')
    area_area = rootgrp.variables['AREA'][:]
    rootgrp.close()
    # delta_mortality in the area of interest
    delta_mort_reg = [
            np.sum(delta_mort[i]*area_area/100) for i in range(pt)]

    # Calculate delta yll in the selected region
    deltayll_reg = [np.sum(delta_yll[i] * area_area/100) for i in range(pt)]
    # yll per 100000 people
    deltayll_spec_reg = [((
            deltayll_reg[i] / np.sum(
                    pop30plus * area_area/100)) * 100000) for i in range(pt)]

    return (deltayll_reg, delta_mort_reg,
            deltayll_spec_reg)

# -------------------------------------------------------------------------
# main program starts here
# -------------------------------------------------------------------------

if __name__ == '__main__':

    from sherpa_globals import (path_tiff, path_mortbaseline,
                                path_dust_conc_cdf_test,
                                path_salt_conc_cdf_test)

    level = 'NUTS_Lv0'
    parea = 'parea'
    path_conc_nc = 'D://sherpa.git//Sherpa//outputdieselgate//delta_concentration.nc'
    codes = ['CY', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES', 'FI',
             'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT', 'NL', 'PL',
             'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
#    codes = ['IT']
    deltayll_reg = {}
    deltayll_spec_reg = {}
    delta_mort_reg = {}
    for code in codes:
        path_areasel_nc = 'workdir\\{}{}.nc'.format(level, code)
        if not os.path.exists(path_areasel_nc):
                area_sel = gridint_toarray(level, parea, code)
                write_nc(area_sel, path_areasel_nc, 'AREA', '%')

        (deltayll_reg[code], delta_mort_reg[code],
         deltayll_spec_reg[code]) = calc_impacts(
                 path_conc_nc, path_areasel_nc, path_tiff, path_mortbaseline,
                 path_dust_conc_cdf_test, path_salt_conc_cdf_test,
                 code, spec='d',
                 approx='l', miller=False, std_life_exp=70)
