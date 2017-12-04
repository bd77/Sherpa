# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017
    
    Health impacts of PM2.5 - WHO-Europe method (HRAPIE reccomendations)
    to calculate mortality according to the Risk Rates
    
    (exact values for beta from the software AirQ+ (WHO))

    Corresonding YLLs or days of life lost are calculated considering the 
    distribution of mortality and population by age and by country, 
    where data is not available in the ICD-10 format YLLs are calculated 
    considering the average value for the countries that 
    are availble.

    NB: A positive delta means a reduction!
    
    - Methodology:
    
    Estimating Local Mortality Burdens associated with Particulate Air 
    Pollution 2014 Public Health England
    
    World Health Organization Europe, 2013. Health risks of air pollution
    in Europe - HRAPIE project - Recommendations for concentration–response
    functions for cost–benefit analysis of particulate matter, ozone and
    nitrogen dioxide, Copenhagen Ø, Denmark.

    Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios
    for the EU Clean Air Package Version. Version 2

    Relative Risk function
    World Health Organization Europe, 2017. AirQ+: software tool for health
    risk assessment of air pollution.
     
    Implementation of the HRAPIE Recommendations for European Air Pollution 
    CBA work 2014 EMRC  (Holland)
      
    Data for baseline population (path_mortbaseline):
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
import json

import os as os

from sherpa_globals import (path_result_cdf_sherpa,
                            path_healthbl, json_path)

def health_impact(pop30plus, pm25_conc, ar_drate, ar_lyl, approx='l'):
    """
    Function that caclulates the health impact
    input : 
        - pop30plus = array with the distribution of the population over 30
        years of age
        - pm25_conc = array with the concentration of PM2.5 (delta or total)
        - ar_drate = array with the distribution of baseline death rate
        (from all cause mortality)
        - ar_lyl = array with the average years of life lost per death 
        over 30 years of age
        - approx = 'e' for exponential and 'l' for linear ('e' is not used)
    output :
        - mort = array with mortality
        - dll = array with the days of life lost per year
        - dll_spec = array with the days of life lost per person per year
    @author: peduzem
    """
# -----------------------------------------------------------------------------
    # create empty arrays to store results
    mort = np.zeros(np.shape(pop30plus))
    dll = np.zeros(np.shape(pop30plus))
    
# -----------------------------------------------------------------------------
    # CONCENTRATION RESPONSE FUNCTION:
    # Estimate of mortality
    # considering bounds for 95% CI
    if approx == 'l':
        # linear approximation
        mrr = 0.06  # 'middle' value
        lrr = 0.04  # lower bound
        hrr = 0.083  # higher bound
        # crf = 0.006 from HRAPIE project
        rr = [lrr, mrr, hrr]
        pt = len(rr)
        mort = mort + (np.where(np.isnan(pm25_conc), 0, (
                     [(rr[i]*pm25_conc/10)/(1+rr[i]*pm25_conc/10)*pop30plus*ar_drate
                     for i in range(len(rr))]))) 
  
    elif approx == 'e':
        # Taken from AirQ+ (WHO)
        # cutoff = 10 # microg/m3 # Taken from AirQ+ (WHO)
        # but not mentioned in the guidelines so not considered here.
        mbeta = 0.006015392281974714  # 'middle' value
        lbeta = 0.003922071315328133  # lower bound
        hbeta = 0.007973496801885352  # higher bound
        beta = [lbeta, mbeta, hbeta]
        pt = len(beta)
        mort = mort + (np.where(
                    np.isnan(pm25_conc), 0, (
                            [(1-(np.exp(-beta[i]*pm25_conc))) *
                             pop30plus * ar_drate
                             for i in range(len(beta))]))) 

# -----------------------------------------------------------------------------
    # ESTIMATE OF THE YLL (Not in the Guidelines)
    # days of life lost per year 
    dll = dll + (np.where(np.isnan(pm25_conc), 0,
                             [mort[i] * ar_lyl * 365
                             for i in range(pt)])) 
    # days of life lost per person per year 
    dll_spec = [np.divide(dll[i], pop30plus, out=np.zeros_like(dll[i]), where=pop30plus!=0) for i in range(pt)] 
    
# -----------------------------------------------------------------------------
    # return results    
    return mort, dll, dll_spec



def main_healthimpact(path_healthbl, path_result_cdf_sherpa, json_path):
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

    # value of concentration model results minus base line
    fh_pm25_natural = Dataset(path_healthbl, mode='r')
    pm25_natural = fh_pm25_natural.variables['conc'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
    fh_pm25_natural.close()
    
    # scenario values minus natural concentration (to check)
    path_value_nc = path_result_cdf_sherpa + 'value_conc.nc'
    fh_pm25_conc = Dataset(path_value_nc, mode='r')
    pm25_conc = fh_pm25_conc.variables['conc'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
    fh_pm25_conc.close()
   
    sce_pm25_conc = pm25_conc - pm25_natural
    
    # delta concentration from model resutls
    path_conc_nc = path_result_cdf_sherpa + 'delta_concentration.nc'
    fh_deltapm25 = Dataset(path_conc_nc, mode='r')
    d_pm25_conc = fh_deltapm25.variables['delta_concentration'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
    fh_deltapm25.close()
    
#    pm25_base = sce_pm25_conc + d_pm25_conc
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
    
    sce_mort, sce_dll, sce_dll_spec = health_impact(pop30plus, sce_pm25_conc, ar_drate, ar_lyl, approx='l')
    delta_mort, delta_dll, delta_dll_spec = health_impact(pop30plus, d_pm25_conc, ar_drate, ar_lyl, approx='l')

              
    dflt_dict = {
    'd_dll': {
            'aggregation': 'sum',
            'ci': ['d_dll_lb', 'd_dll', 'd_dll_ub'],
            'combo_box': 'delta days of life loss',
            'long_description': ['delta days of life loss lower bound',
                                 'delta days of life loss',
                                 'delta days of life loss upper bound'],
            'units': 'dll/year'},
    'd_dll_pp': {
            'aggregation': 'population weighted average',
            'ci': ['d_dll_pp_lb', 'd_dll_pp', 'd_dll_pp_ub'],
            'combo_box': 'delta days of life loss',
            'long_description': 
                ['delta days of life loss per person lower bound',
                 'delta days of life loss per person',
                 'delta days of life loss per person upper bound'],
            'units': 'dll/(person year)'},
    'd_mort': {
            'aggregation': 'sum',
            'ci': ['d_mort_lb', 'd_mort', 'd_mort_ub'],
            'combo_box': 'delta mortality',
            'long_description': ['delta mortality lower bound',
                                 'delta mortality',
                                 'delta mortality upper bound'],
            'units': 'people/year'},
    'v_dll': {
            'aggregation': 'sum',
            'ci': ['v_dll_lb', 'v_dll', 'v_dll_ub'],
            'combo_box': 'days of life loss',
            'long_description': ['days of life loss lower bound',
                                 'days of life loss',
                                 'days of life loss upper bound'],
            'units': 'dll/year'},
    'v_dll_pp': {
            'aggregation': 'population weighted average',
            'ci': ['v_dll_pp_lb', 'v_dll_pp', 'v_dll_pp_ub'],
            'combo_box': 'days of life loss per person',
            'long_description': ['days of life loss per person lower bound',
                                 'days of life loss per person',
                                 'days of life loss per person upper bound'],
            'units': 'dll/(person year)'},
  'v_mort': {
          'aggregation': 'sum',
          'ci': ['v_mort_lb', 'v_mort', 'v_mort_ub'],
          'combo_box': 'mortality',
          'long_description': ['mortality lower bound',
                               'mortality',
                               'mortality upper bound'],
          'units': 'people/year'}}

 
    json_path = 'config\\config.json'
    if os.path.exists(json_path):    
        print('Using stored json file')
        json_file = open(json_path)
        json_str = json_file.read()
        cfg_dct = json.loads(json_str)
    else:
        cfg_dct = dflt_dict
        
    outfile=path_result_cdf_sherpa + 'healthimp.nc'
    if os.path.exists(outfile):
        os.remove(outfile)   
    for key in cfg_dct.keys():
        if key == 'd_mort':
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(delta_mort[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_mort': 
            for it in enumerate(cfg_dct[key]['ci']):
                write_nc(sce_mort[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'd_dll':
            for it in enumerate(cfg_dct[key]['ci']):
                write_nc(delta_dll[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_dll': 
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(sce_dll[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'd_dll_pp':
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(delta_dll_spec[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_dll_pp': 
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(sce_dll_spec[it[0]], outfile, it[1], cfg_dct[key]['units'],
                     addnutsid=True, l_name=cfg_dct[key]['long_description'][it[0]])
    
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
    rootgrp = Dataset(path_healthbl, 'r')
    lon_array = rootgrp.variables['longitude'][:]
    lat_array = rootgrp.variables['latitude'][:]
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
    
    main_healthimpact(path_healthbl, path_result_cdf_sherpa, json_path)


