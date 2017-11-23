# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017
    
    Health impacts of PM2.5 - WHO-Europe method (HRAPIE reccomendations)
    to calculate mortality according to the Concentration-Response-Function
    (exact values from the software AirQ+ (WHO))

    Corresonding YLLs or days of life lost are calculated considering the 
    distribution of mortality and population by age and by country, 
    where data is not available in the ICD-10 format YLLs are calculated 
    considering the average value for the countries that 
    are availble.

    NB: A positive delta means a reduction!
    
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


    @author: peduzem
    """


from netCDF4 import Dataset  
import numpy as np

import os as os

from sherpa_globals import (path_base_conc_cdf_test, path_dust_conc_cdf_test,
                            path_salt_conc_cdf_test, path_result_cdf_test,
                            path_healthbl, path_model_cdf_test)

def health_impact(pop30plus, pm25_conc, ar_drate, ar_lyl, approx='l'):
    """
    Function that caclulates the health impact
    input : 
        - pop30plus = array with the distribution of the population over 30
        years of age
        - pm25_conc = array with the concentration of PM2.5
        - ar_drate = array with the distribution of baseline death rate
        (from all cause mortality)
        - ar_lyl = array with the average years of life lost per death 
        over 30 years of age
        - approx = 'e' for exponential and 'l' for linear
    output :
        - delta_mort = array with mortality
        - delta_dll = array with the days of life lost per year
        - detla_dll_spec = array with the days of life lost per person per year
    @author: peduzem
    """
# -----------------------------------------------------------------------------
    # create empty arrays to store results
    delta_mort = np.zeros(np.shape(pop30plus))
    delta_dll = np.zeros(np.shape(pop30plus))
    
# -----------------------------------------------------------------------------
    # CONCENTRATION RESPONSE FUNCTION:
    # Estimate of mortality
    # considering bounds for 95% CI
    if approx == 'l':
        # linear approximation
        mrr = 1.06  # 'middle' value
        lrr = 1.04  # lower bound
        hrr = 1.083  # higher bound
        # crf = 0.006 from HRAPIE project
        af = np.asarray([(lrr-1)/lrr, (mrr-1)/mrr, (hrr-1)/hrr]) / 10
        pt = len(af)
        delta_mort = delta_mort + (np.where(np.isnan(pm25_conc), 0, (
                     [af[i]*pm25_conc*pop30plus*ar_drate
                     for i in range(len(af))]))) 
    elif approx == 'e':
        # Taken from AirQ+ (WHO)
        # cutoff = 10 # microg/m3 # Taken from AirQ+ (WHO)
        # but not mentioned in the guidelines so not considered here.
        mbeta = 0.006015392281974714  # 'middle' value
        lbeta = 0.003922071315328133  # lower bound
        hbeta = 0.007973496801885352  # higher bound
        beta = [lbeta, mbeta, hbeta]
        pt = len(beta)
        delta_mort = delta_mort + (np.where(
                    np.isnan(pm25_conc), 0, (
                            [(1-(np.exp(-beta[i]*pm25_conc))) *
                             pop30plus * ar_drate
                             for i in range(len(beta))]))) 

# -----------------------------------------------------------------------------
    # ESTIMATE OF THE YLL (Not in the Guidelines)
    # days of life lost per year 
    delta_dll = delta_dll + (np.where(np.isnan(pm25_conc), 0,
                             [delta_mort[i] * ar_lyl * 365
                             for i in range(pt)])) 
    # days of life lost per person per year 
    delta_dll_spec = [np.divide(delta_dll[i], pop30plus, out=np.zeros_like(delta_dll[i]), where=pop30plus!=0) for i in range(pt)] 
    
# -----------------------------------------------------------------------------
    # return results    
    return delta_mort, delta_dll, delta_dll_spec



def main_healthimpact(path_base_conc_cdf_test, path_dust_conc_cdf_test, path_salt_conc_cdf_test, path_healthbl, path_result_cdf_test):
    """
    Main functin that calculates the health impacts given the paths: 
    input: 
        - path_base_conc_cdf_test = base case concentration 
        - path_dust_conc_cdf_test = path of baseline dust concentration 
        - path_salt_conc_cdf_test = path of baseline salt concentration 
        - path_healthbl = path where results are stored (health baseline)
        - path_result_cdf_test: path of the delta concentrations
           (output of module1) 
    """
    # base case PM25 conentration
    rootgrp = Dataset(path_base_conc_cdf_test, mode='r')
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

    # delta concentration from model resutls
    path_conc_nc = path_result_cdf_test + 'delta_concentration.nc'
    fh_deltapm25 = Dataset(path_conc_nc, mode='r')
    d_pm25_conc = fh_deltapm25.variables['delta_concentration'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
    fh_deltapm25.close()
    
    # get baseline data from nc file
    fh = Dataset(path_healthbl, mode='r')
    pop30plus = fh.variables['ppl30+'][:]
    fh.close()
    fh = Dataset(path_healthbl, mode='r')
    ar_drate = fh.variables['deathsppl30+'][:]
    fh.close()
    fh = Dataset(path_healthbl, mode='r')
    ar_lyl = fh.variables['lyl30+'][:]
    fh.close()
    
    # calculate impacts
    delta_mort, delta_dll, delta_dll_spec = health_impact(pop30plus, pm25_conc, ar_drate, ar_lyl, approx='l')

    # file to write resutls (remove if it already exists)
    outfile=path_result_cdf_test + 'healthimp.nc'
    if os.path.exists(outfile):
        os.remove(outfile)   
    
    # write all results into netcdf files
    write_nc(delta_mort[1], outfile, 'bl_mort', 'number', addnutsid=True, l_name='base line mortality')
    write_nc(delta_mort[0], outfile, 'bl_mort_lb', 'number', addnutsid=True, l_name='base line mortality lower bound')
    write_nc(delta_mort[2], outfile, 'bl_mort_ub', 'number', addnutsid=True, l_name='base line mortality upper bound')
    write_nc(delta_dll[1], outfile, 'bl_dll', 'dll per year', addnutsid=True, l_name='base line days of life loss per year')
    write_nc(delta_dll[0], outfile, 'bl_dll_lb', 'dll per year', addnutsid=True, l_name='base line days of life per year loss lower bound')
    write_nc(delta_dll[2], outfile, 'bl_dll_ub', 'dll per year', addnutsid=True, l_name='base line days of life per year loss upper bound')
    write_nc(delta_dll_spec[1], outfile, 'bl_dll_pp', 'dll per person per year', addnutsid=True, l_name='base line days of life loss per year')
    write_nc(delta_dll_spec[0], outfile, 'bl_dll_pp_lb', 'dll per person per year', addnutsid=True, l_name='base line days of life loss per person per year lower bound')
    write_nc(delta_dll_spec[2], outfile, 'bl_dll_pp_up', 'dll per person per year', addnutsid=True, l_name='base line days of life loss per person per year upper bound')
    
    delta_mort, delta_dll, delta_dll_spec = health_impact(pop30plus, d_pm25_conc, ar_drate, ar_lyl, approx='l')
    write_nc(delta_mort[1], outfile, 'mort', 'number', addnutsid=True, l_name='mortality')
    write_nc(delta_mort[0], outfile, 'mort_lb', 'number', addnutsid=True, l_name='mortality lower bound')
    write_nc(delta_mort[2], outfile, 'mort_ub', 'number', addnutsid=True, l_name='mortality upper bound')
    write_nc(delta_dll[1], outfile, 'dll', 'dll per year', addnutsid=True, l_name='days of life loss per year')
    write_nc(delta_dll[0], outfile, 'dll_lb', 'dll per year', addnutsid=True, l_name='days of life loss per year upper bound')
    write_nc(delta_dll[2], outfile, 'dll_up', 'dll per year', addnutsid=True, l_name='days of life loss per year')
    write_nc(delta_dll_spec[1], outfile, 'dll_pp', 'dll per person per year', addnutsid=True, l_name='days of life loss per person per year')
    write_nc(delta_dll_spec[0], outfile, 'dll_pp_lb', 'dll per person per year', addnutsid=True, l_name='days of life loss per person per year lower bound')
    write_nc(delta_dll_spec[2], outfile, 'dll_pp_up', 'dll per person per year', addnutsid=True, l_name='days of life loss per person per year upper bound')

## SUPPORT FUNCTIONS (IDEALLY IN THE AUXIALIARIES FILE)

def write_nc(array, path_nc, name_var, unit_var, addnutsid=False, l_name=None):
    ''' Function to write an array in a netcdf file,
        input:
            - array: data to write
            - path_nc: path of netcdf file
            - name_var: name for data in array
            - unit_var: units for data in array
            - addnutsid: if True the layer nuts_id is added so that the
                nectcdf file is consistent with the ones provided
                by terraria
    @author: peduzem
    '''
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()
    if not os.path.exists(path_nc):
        mode = 'w' 
        fh=Dataset(path_nc, mode=mode, format='NETCDF3_CLASSIC') 
        fh.createDimension('latitude', len(lat_array))
        fh.createDimension('longitude', len(lon_array))
        latitude = fh.createVariable('latitude', 'f4', ('latitude',))
        longitude = fh.createVariable('longitude', 'f4', ('longitude',)) 
        if addnutsid is True:
#        fh.createDimension('z', 10)
            fh.createDimension('nuts_id', 1)
            var = fh.createVariable(name_var, 'f8',
                                    ('nuts_id', 'latitude', 'longitude',))
            nutsid = fh.createVariable('NUTS', 'i4', ('nuts_id',))
            longitude[:] = lon_array
            latitude[:] = lat_array
            nutsid[0] = 1
            var[0, :] = array
        elif addnutsid is False:
            longitude[:] = lon_array
            latitude[:] = lat_array
            var = fh.createVariable(name_var, 'f8', ('latitude', 'longitude'))
            var[:] = array          
    else:
        mode = 'a'
        fh=Dataset(path_nc, mode=mode, format='NETCDF3_CLASSIC')
        if addnutsid is True:
            var = fh.createVariable(name_var, 'f8',
                                    ('nuts_id', 'latitude', 'longitude',))
            var[0, :] = array
        elif addnutsid is False:
            var = fh.createVariable(name_var, 'f8', ('latitude', 'longitude'))
            var[:] = array

    fh.variables[name_var].units = unit_var
    if l_name is not None:
        fh.variables[name_var].long_name =l_name
    fh.close()  

if __name__ == '__main__':
    
    main_healthimpact(path_base_conc_cdf_test, path_dust_conc_cdf_test, path_salt_conc_cdf_test, path_healthbl, path_result_cdf_test)
    
    level = 'NUTS_Lv0'
    code = 'AT'
    path_areasel_nc = 'workdir\\{}{}f8.nc'.format(level, code)
    path_healthres = path_result_cdf_test + 'healthimp.nc'
    delta_mort = {}
    rootgrp = Dataset(path_healthres, 'r')
    delta_mort[0] = rootgrp.variables['bl_mort_lb'][:]
    delta_mort[1] = rootgrp.variables['bl_mort'][:]
    delta_mort[2] = rootgrp.variables['bl_mort_ub'][:]
    rootgrp.close()
    
    rootgrp = Dataset(path_areasel_nc, mode='r')
    area_area = rootgrp.variables['AREA'][:]
    rootgrp.close()
    # delta_mortality in the area of interest
    pt = 3
    delta_mort_reg = [
            np.sum(delta_mort[i]*area_area/100) for i in range(pt)]
    area = np.sum(area_area)
    area2 = np.sum(area_co['AT'])