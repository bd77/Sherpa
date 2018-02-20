# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017
    This script generates the input file for the health impact assessment:
    moduel8_healthia.py
        
    This concerns only the health impacts of PM2.5 - WHO-Europe method 
    (HRAPIE reccomendations)

    MAIN ASSUMPTIONS:
        
        - BASELINE VALUES ARE HOMOGENEOUS WITHIN A COUNTRY AND AVERAGED BETWEEN
          COUNTRIES FOR BORDER CELLS!! 
        - RISK RATE FUNCTIONS ARE APPLIED ONLY TO ANTHROPOGENIC PM2.5 as in
          EC4MACS Modelling Methodology the ALPHA Benefit Assessment  Model and 
          Service Contract for Carrying out Cost-Benefit Analysis of Air Quality 
          Related Issues, in particular in the Clean Air for Europe (CAFE)
          Programme.
          In CAFE: 
          " The case is sometimes made that there are natural backgrounds of 
          ozone (and of other pollutants also) and that either:        
              - i. There are no adverse health effects associated with 
                concentrations below these backgrounds, because they are natural
              - ii. Any associated adverse health effects should not be quantified,
                because it is impossible to reduce pollution to below these levels"

    Required files:
        - path_mortbaseline: excel file with the data of the baseline
        - path_tiff: .tif file with the population distribution (LUISA)
          (For EU28 data is from LUISA, in its last version 11/22/2017,
          relative to the year 2015. Extra EU, always for 2015, is taken from:
          http://sedac.ciesin.columbia.edu/data/collection/gpw-v4
        - grid intersect files 
          Country codes for the grid intesect are:
             'CY', 'ES', 'HU', 'ME', 'NL', 'SI', 'AT', 'BE', 'BG', 'CH', 
             'CZ', 'DE', 'DK', 'EE', 'EL', 'FI', 'FR', 'HR', 'IE', 'IT',
             'LI', 'LT', 'LU', 'LV', 'MK', 'MT', 'NO', 'PL', 'PT', 'RO', 
             'SE', 'SK', 'TR', 'UK'.
        - path_model_cdf_test: (only to be able to save ncdf with the same 
          format)

    BIBLIOGRAPHY:

    - Methodology:
        World Health Organization Europe, 2013. Health risks of air pollution
        in Europe - HRAPIE project - Recommendations for concentration–response
        functions for cost–benefit analysis of particulate matter, ozone and
        nitrogen dioxide, Copenhagen Ø, Denmark.

        Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios
        for the EU Clean Air Package Version. Version 2

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


from netCDF4 import Dataset  # for using netCDF files
import numpy as np
# from osgeo import gdal, ogr, osr  # conda install -c conda-forge gdal
import pandas as pd
# import scipy.io as sio
import os as os
from osgeo import gdal
import sys
import time # time the code, check which are the bottle necks

from sherpa_globals import (path_tiff, path_mortbaseline, path_grid_txt,
                            path_healthbl, path_model_cdf_test, 
                            path_dust_conc_cdf_test, 
                            path_salt_conc_cdf_test)

def baseline_nc(path_tiff, path_mortbaseline, path_healthbl, 
                path_dust_conc_cdf_test, 
                path_salt_conc_cdf_test, path_model_cdf_test, std_life_exp=70):
    """
    Function that produces the base line netcdf for the SHERPA tool 
    
    input : 
        - path_tiff = path to the tiff image with the population from Marco
        - path_mortbaseline = path to the mortality baseline file with 
          the information per country
        - path_dust_conc_cdf_test = path of baseline dust concentration 
        - path_salt_conc_cdf_test = path of baseline salt concentration 
        - path_healthbl = path where results are stored (health baseline)
        - std_life_exp = standard life expectancy to calculate the years of
          life loss, 70 is the value used by WHO in
          http://data.euro.who.int/dmdb/
          
    output :
        - healthbl_nc.nc' = netcdf file with: 
            - ppl30+ = population over 30 distribution
            - deathsppl30+ = death rate for people over 30 distribution 
            - lyl30+ = life of years lost for each person dying over 30. 
            - pm25_natural = natural PM (from dust and salt)
            
    @author: peduzem
    """
# -----------------------------------------------------------------------------
    # remove output file is one already exists
    if os.path.exists(path_healthbl):
        os.remove(path_healthbl)
# -----------------------------------------------------------------------------
    # BASELINE POPULATION DATA and other PRECOMPUTATION
    # Load population file from LUISA
    popall = tiftogridgeneral(path_tiff)
    # remove all negative and -inf values (not sure why they are there)
    popall = np.where(np.isinf(popall), 0, popall)
    popall = np.where(popall < 0, 0, popall)
#    path_pop = 'D://sherpa.git//Sherpa//input//pop.nc'
#    write_nc(popall, path_pop, 'pop', 'number', path_model_cdf_test, addnutsid=True, l_name='population')

#    popall = np.where(np.negative(popall), 0, popall)
    popall.sum()
    # READ EXCEL FILE WITH BASELINE DATA FROM THE WHO
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
        because the CRF does not depend on age.
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
    
    # Calculting average values to use if the country selected is not in
    # the database. Average values calculating considering EU28 countries
    # which are also in the WHO database

    countries = ['AT', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL',  
                 'ES', 'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LI', 'LT', 'LU', 
                 'LV', 'ME', 'MK', 'MT', 'NL', 'NO', 'PL', 'PT', 'RO', 'SE',
                 'SI', 'SK', 'TR', 'UK']
    
    # Find coutnries that are both in the 34 countries considered in SHERPA 
    # and in the mortality database
    intersect = list(
            set.intersection(set(countries), set(list(df_mort.index))))
   
    # Calculate averages for the countries in both lists:
    # these values are used when WHO data is not available
    m_drate = drate.ix[intersect].mean(axis=0)
    m_p30_ptot = p30_ptot.ix[intersect].mean(axis=0)
    m_lyl = df_lyl.ix[intersect].mean(axis=0)

    # Calculate the percentage of each cell belonging to each country
    level = 'NUTS_Lv0'
    parea = 'parea'
    area_co = {}
    ttot = 0
    for country in countries:
#        country = 'AT'
        t0 = time.time()
        area_co[country] = gridint_toarray(level, parea, country, path_model_cdf_test)
        t1 = time.time()
        dt = t1 - t0
        print('to calculate area of country', country, 'it takes', dt)
        ttot = ttot + dt
     
    # preparing arrays to store results:
    pop30plus = np.zeros(np.shape(popall))
    ar_drate = np.zeros(np.shape(popall))
    ar_lyl = np.zeros(np.shape(popall))
#    ar_le = np.zeros(np.shape(popall))
    
    for country in countries:
        if country in list(df_mort.index):
            # Scaling factor to respect country totals as provided by WHO
            # in each grid cell of the area of interest
            # warning: this is valid only in the country
            pop30plus = pop30plus + (np.array(popall * p30_ptot.ix[country])
                                     * area_co[country] / 100)
            ar_drate = ar_drate + drate.ix[country] * area_co[country] / 100
            ar_lyl =  ar_lyl + df_lyl.ix[country] * area_co[country] / 100  
#            ar_le  = ar_le + df_le.ix[country][2009] * area_co[country] / 100                           
        else:
            # Scaling to only the population over 30 in each grid cell
            # warning: this is valid only in the country
            pop30plus = pop30plus + (np.array(popall * m_p30_ptot) 
                        * area_co[country] / 100)            
            ar_drate = ar_drate + m_drate * area_co[country] / 100
            ar_lyl =  ar_lyl +  m_lyl * area_co[country] / 100
#            ar_le  = ar_le + m_le * area_co[country] / 100  

                                            
                                            
#    # base case PM25 conentration
#    rootgrp = Dataset(path_base_conc_cdf_test, mode='r')
#    bc_pm25_conc = rootgrp.variables['conc'][:]
##    bc_pm25_units = rootgrp.variables['conc'].units
#    rootgrp.close()

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
    pm25_natural = pSALT25 + pDUST25
    
    
    
    # write resutls in nc files
    write_nc(pm25_natural, path_healthbl, 'conc', u'\u03BC'+'g/m'+ u'\u00B3', path_model_cdf_test, addnutsid=True, l_name='natural (salt and dust) PM2.5 concentration')
    write_nc(pop30plus, path_healthbl, 'ppl30+', '-', path_model_cdf_test, addnutsid=True, l_name='number of people avove 30 years of age')
    write_nc(ar_drate, path_healthbl, 'deathsppl30+', '-', path_model_cdf_test, addnutsid=True, l_name='mortality rate - people above 30 over number of death above 30')
    write_nc(ar_lyl, path_healthbl, 'lyl30+', '-', path_model_cdf_test, addnutsid=True, l_name='average years of life lost for each death above 30')           
# -----------------------------------------------------------------------------



## SUPPORT FUNCTIONS (IDEALLY IN THE AUXIALIARIES FILE)

def write_nc(array, path_nc, name_var, unit_var, path_model_cdf_test, addnutsid=False, l_name=None):
    ''' Function to write an array in a netcdf file,
        input:
            - array: data to write
            - path_nc: path of netcdf file
            - name_var: name for data in array
            - unit_var: units for data in array
            - path_model_cdf_test: path of model to copy lat and lon array from 
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
            var = fh.createVariable(name_var, 'f4',
                                    ('nuts_id', 'latitude', 'longitude',))
            nutsid = fh.createVariable('NUTS', 'i4', ('nuts_id',))
            longitude[:] = lon_array
            latitude[:] = lat_array
            nutsid[0] = 1
            var[0, :] = array
        elif addnutsid is False:
            longitude[:] = lon_array
            latitude[:] = lat_array
            var = fh.createVariable(name_var, 'f4', ('latitude', 'longitude'))
            var[:] = array          
    else:
        mode = 'a'
        fh=Dataset(path_nc, mode=mode, format='NETCDF3_CLASSIC')
        if addnutsid is True:
            var = fh.createVariable(name_var, 'f4',
                                    ('nuts_id', 'latitude', 'longitude',))
            var[0, :] = array
        elif addnutsid is False:
            var = fh.createVariable(name_var, 'f4', ('latitude', 'longitude'))
            var[:] = array

    fh.variables[name_var].units = unit_var
    if l_name is not None:
        fh.variables[name_var].long_name =l_name
    fh.close()

def tiftogridgeneral(path_tiff):
    gdal.UseExceptions()
    ds = None
    try:
        ds = gdal.Open(path_tiff)
        em = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        # re-arrange emission inventory file for Sherpa's grid
        # initialize array
        Amine = em  # Amine : matrix
        for i in range(0, 382):
            ind1 = 2*(i-1)  # included
            ind2 = 1+2*(i-1)+1  # excluded
            Amine[:, i-1] = (np.sum(em[:, ind1:ind2], axis=1))
        Amine[:, 382:384] = 0
        # Cancelling the extra columns and extra rows
        # (is there a better way to do this?**)
        for deli in range (0,144):
            Amine = np.delete(Amine, (0), axis=0) # delete first 144 rows
        for delj in range (0,398): # delete last 398 columns
            Amine = np.delete(Amine, (383), axis=1)
        Amine_T = Amine[np.newaxis]
        Afinal=np.fliplr(Amine_T)
        arr=Afinal
        return arr
    except(RuntimeError, AttributeError):
        pass
    
#def read_nc(nc_file):
#    '''
#    NAME
#        Reads SHERPA ncdf file with Python
#    PURPOSE
#        To read matrix data and put them in a multindexed dataframe
#    PROGRAMMER(S)
#        Denise Pernigotti
#    REVISION HISTORY
#    
#    REFERENCES
#    
#    '''
#    #nc_file='input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
#    #nc_file='input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
##    nc_file = 'D:/programs/sherpa/app/data/temp/delta_emission.nc'
#    nc_data = Dataset(nc_file, 'r')
#    nc_dims = [dim for dim in nc_data.dimensions]
#    nc_vars = [var for var in nc_data.variables]
#    #sometimes the latitude is written just with lat as in model data
#    latname=list(filter(lambda x: x in nc_vars, ['Lat','lat','latitude']))[0]
#    #latname=list(filter (lambda x: 'lat' in x, nc_vars))[0]
#    lats = nc_data.variables[latname][:]
#    lonname=list(filter(lambda x: x in nc_vars, ['Lon','lon','longitude']))[0]
#    #lonname=list(filter (lambda x: 'lon' in x, nc_vars))[0]
#    lons = nc_data.variables[lonname][:]
#    #if there are three dimensional arrays
#    if len(nc_dims)==3:
#        ncz=str(list(set(nc_dims)-set(['latitude','longitude']))[0])
#        nz=range(len(nc_data.dimensions[ncz]))
#        if ncz=='pollutant':
#            strpoll=nc_data.Order_Pollutant
#            nznames=strpoll.split(', ')
#        else:
#            nznames=[ncz +"{:02d}".format(x+1) for x in nz]
#            #nznames=[ncz + s for s in map(str,range(1,len(nc_data.dimensions[ncz])+1))]
#    #create an index with lat and lon
#    #latrep=map(str, np.repeat(lats,len(lons)))
#    #lonrep=map(str, np.tile(lons,len(lats)))
#    #trasform variables arrays in vectors
#    #allvar={'lat_lon':map(lambda (x,y): x+'_'+y, zip(latrep, lonrep))}
#    #create lat and lon info
#    if len(lats.shape)==2 and len(lons.shape)==2:
#        nrow=lats.shape[0]
#        ncol=lats.shape[1]
#        lon=lons.ravel()
#        lat=lats.ravel()
#    else:
#        nrow=len(lats)
#        ncol=len(lons)
#        lon=np.tile(lons,nrow)
#        lat=np.repeat(lats,ncol)
#
#    y=np.repeat(range(1, nrow+1),ncol)
#    x=np.tile(range(1, ncol+1),nrow)
#    row=list(map(str,y))
#    col=list(map(str,x))
#    index_grid=list(map(lambda x: '_'.join(x),list(zip(col,row))))
#
#    allvar={}
#    allvar['coord']=pd.DataFrame(lon,columns=['lon'])
#    allvar['coord']['lat']=lat
#    allvar['coord']['x']=x
#    allvar['coord']['y']=y
#    allvar['coord'].index=index_grid
#    nc_vars.remove(latname)
#    nc_vars.remove(lonname)
#    # added by EPE to deal with the delta_emissions file created byt 
#    # the GUI which has an extra variable 
#    # @todo this condition can be removed in the future if the GUI 
#    # is fixed (20180219)
#    if 'Nsnaps' in nc_vars: 
#        nc_vars.remove('Nsnaps')
#        
#    for var in nc_vars:
#        varnc=nc_data.variables[var][:]
#        if len(nc_dims)==3:
#            allvar[var]=pd.concat(map(lambda sn : pd.Series(varnc[sn].ravel()),nz),axis=1)
#            allvar[var].columns=nznames
#        else:
#            allvar[var]=pd.DataFrame(varnc.ravel())
#            allvar[var].columns=[var]
#        allvar[var].index=index_grid
#        #allvarnc[var]=allvarnc[var].transpose()
#        #index_var = pd.MultiIndex.from_tuples(zip(np.repeat(var,len(nz)),nznames), names=['vars', ncz])
#        #allvar[var].columns=index_var
#    nc_data.close()
#    reform = {(outerKey, innerKey): values for outerKey, innerDict in allvar.items() for innerKey, values in innerDict.items()}
#    df=pd.DataFrame(reform)
#    return df.transpose()

def read_nuts_area(filenuts,calcall=False,nullnut=None,nutsall=None):
    '''
    NAME
        Import info on grid points attribution to nuts or specific area type from ascii file
    PURPOSE
        Import info on grid points attribution to nuts or specific area type from ascii file/s.
        If the file is single then it must contain the column 'Area [km2]' relative to % of the area in the finest nut,
        this datum will be set to each nut but it will then aggregated for larger nuts when nutsarea will be calculated
        If the files are two, then each nut will have its own % area for each grid point, then the data will be merged here
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
    
    REFERENCES
    
    '''
    nuts_info_all={}
    if(filenuts != 'rect'):
        nuts_def= filenuts +'.txt'
        nuts_info = pd.read_csv(nuts_def,delimiter="\t")
        nuts_info=nuts_info.dropna(axis=1,how='all')
        nutsnames=list(nuts_info.columns[~nuts_info.columns.isin(['POP','COL','ROW','Area [km2]','LAT','LON'])])
        #optional 'nut' comprising all grid points
        if calcall :
        #nutsnames.insert(0, 'ALL')
            nutsnames.insert(0, 'ALL_'+nutsnames[0])
            nuts_info[nutsnames[0]]=nutsnames[0]
        nuts_info['grid']=['_'.join(str(i) for i in z) for z in zip(nuts_info['COL'],nuts_info['ROW'])]
        if 'Area [km2]' in nuts_info.columns:
            nuts_area=pd.concat(map(lambda p: nuts_info['Area [km2]'],nutsnames),axis=1)
            #nuts_area.index=nuts_info['grid']
            nuts_area.columns=nutsnames
           #nuts_info=nuts_info[nutsnames]
        else:
            sys.exit("missing infos on grid cells area per nut")

        #aggregate data for each nut, create a dictionary
        for nut in nutsnames:
            #create a multindex
            index = pd.MultiIndex.from_tuples(list(zip(nuts_info[nut],nuts_info['grid'])), names=['nutname','grid'])
            nut_info=pd.Series(list(nuts_area[nut]), index=index)
            nut_info=nut_info.to_frame(name='area')
            #aggregate data on these nuts if necessary
            nut_info_nut=nut_info.groupby(level=[0,1]).sum()
            #find total area
            grid_area_tot=nut_info_nut.groupby(level=['grid']).sum()
            nut_info_nut['parea']=nut_info_nut/grid_area_tot
            nut_info_nut.loc[nut_info_nut['area']==0,'parea']=0.
            #eventually remove the fillng code
            if nullnut is not None:
                nut_info_nut=nut_info_nut.drop(nullnut, level='nutname')
            nuts_info_all[nut]=nut_info_nut
    else:
        nuts_rect=nutsall
        nuts_rect.index=nuts_rect.index.droplevel(level=0)
        grid_inrect=nuts_rect.index
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']>=rect_coord['ll']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']<=rect_coord['ur']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']>=rect_coord['ll']['lat']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']<=rect_coord['ur']['lat']]
        nuts_rect=nuts_rect.loc[list(grid_inrect)]
        nuts_rect['nutname'] = 'rect'
        nuts_rect.set_index('nutname', append=True, inplace=True)
        nuts_info_all['rect']=nuts_rect.swaplevel(i=-2, j=-1, axis=0)
    return nuts_info_all



def gridint_toarray(level, parea, code, path_model_cdf_test):
    """
    Reads the grid intersect txt files and creates an array with the specified
    dimensions with the percentage of each cell beleonging to the specified area
  . Needs to import the function from DENISEP
    from module7_custom import read_nuts_area

    INPUT:
        level = 'NUTS_Lv0' # or GCITY_CODE or FUA_CODE
        parea = 'parea'
        code = 'IT' or corrisponding codes

    OUTPUT:
        area_sel : array of the percentage (0-100) of each cell belonging
        to the selected area

    @author: peduzem
    """
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()

    nuts_info = read_nuts_area(path_grid_txt, calcall=True)
#    nuts_info.update(read_nuts_area(gcities_txt, nullnut='LAND000'))
#    nuts_info.update(read_nuts_area(fua_txt, nullnut='LAND000'))

    # get area of interest
    narea = nuts_info[level][parea][code]

    # get rows and columns indices
    cols = [list(narea.keys())[i].split("_", 1)[1]
            for i in range(0, len(list(narea.keys())))]
    rows = [list(narea.keys())[i].split("_", 1)[0]
            for i in range(0, len(list(narea.keys())))]

    # convert them from str to int
    rowsint = [int(i) for i in rows]
    colsint = [int(i) for i in cols]

    # create a matrix of the same size as sherpas grid
    area_sel = np.zeros((1, len(lat_array), len(lon_array)))

    # get the points that are not zeros from the gridintersect, columns are
    # the x and rows are the y!
    points = list(zip(np.zeros(len(colsint), dtype=int), colsint, rowsint))

    # assign the values of area fraction (parea)
    # roll both axis by one as the indeces in the grid interesect start from
    # one
    for point in points:
        area_sel[point[0], (point[1]-1), (point[2]-1)] = (
                 narea.get('%d_%d' % ((int(point[(2)])), int(point[(1)]))))

    return area_sel*100

if __name__ == '__main__':

    baseline_nc(path_tiff, path_mortbaseline, path_healthbl, 
                path_dust_conc_cdf_test, 
                path_salt_conc_cdf_test, path_model_cdf_test, std_life_exp=70)
    
#    main_healthimpact(path_base_conc_cdf_test, path_dust_conc_cdf_test, path_salt_conc_cdf_test, path_healthbl, path_result_cdf_test)

