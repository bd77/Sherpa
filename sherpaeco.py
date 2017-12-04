# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:26:42 2016
@author: peduzem

Module to calculate the impact of measures

INPUT:
     *to be completed/updated*
     Need to add an input directory with the following files:
        input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc: using it to get the
        values of lat and lon, will not be necessary in the future
        input/population.mat: population file by Marco
        # input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc' : base case conc
        input/EMI_RED_ATLAS_NUTS0.nc: or NUTS1 to select the area of interest
        input/JRC01.nc : area of each cell
        input/ef_reduction_sherpaeco : csv file with the reduction per
        sector/fuel (see example)

OUTPUT:
    *to be completed/updated*
    input/sherpaeco_reduction.txt: file of emission reductions per macro sector
    Delta Years of Life Lost and Delta Mortality in the area of interest
    and elsewhere
    deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg
    as for module1:
    netcdf with concentration changes per pollutant and cell
    delta emission netcdf with emission changes per precursor and cell
... (to be continued)
"""
import numpy as np
import pandas as pd
import os as os
import os.path

from osgeo import gdal  # ogr, osr
from netCDF4 import Dataset
from os import remove  # for module1
from time import time  # for module1

# for plotting
# from mpl_toolkits.basemap import Basemap
# import matplotlib.pyplot as plt
# import plotly.plotly as py
# from PIL import Image
# import scipy.io as sio

from sherpa_auxiliaries import (create_emission_reduction_dict,
                                create_emission_dict, create_window,
                                read_progress_log, write_progress_log)
from sherpa_auxiliaries_epe import (save_obj, load_obj, gridint_toarray,
                                    write_nc)

from sherpa_globals import (path_area_cdf_test, path_model_cdf_test,
                            path_emission_cdf_test, path_result_cdf_test,
                            path_base_conc_cdf_test, path_mortbaseline, 
                            path_dust_conc_cdf_test, path_salt_conc_cdf_test, 
                            path_tiff)

import module1 as shrp
import module1ema as shrp1
from healthia import calc_impacts


#  @TODO upadate with dataframe function in auxiliaries epe
def write_reductions(path_reduction_txt, red):
    """
    Function to write the reduction file, that is the percentage reduction
    per pollutant per macrosector.
    At the moment only for MS7 (**will be extended to the other sectors)

    Input:
        - path_reduction_txt: path to the reduction file
          (** can be from sherpa_globals)
        - red: array with the reduction (>0) % per pollutant
          (** the order is important)
    Output:
        - path_retduction_txt.txt: with the % reduction of pollutants of MS7

    @author: peduzem
    """
    text_file = open(path_reduction_txt, "w")
    text_file.write("POLL	MS1	MS2	MS3	MS4	MS5	MS6	MS7	MS8	MS9	MS10\n")
    text_file.write("NOx	0	0	0	0	0	0	%s	0	0	0\n" % red[0])
    text_file.write("NMVOC	0	0	0	0	0	0	%s	0	0	0\n" % red[1])
    text_file.write("NH3	0	0	0	0	0	0	%s	0	0	0\n" % red[2])
    text_file.write("PPM	0	0	0	0	0	0	%s	0	0	0\n" % red[3])
    text_file.write("SOx	0	0	0	0	0	0	%s	0	0	0\n" % red[4])
    text_file.close()


# -----------------------------------------------------------------------------
#  @TODO upadate with general function in auxiliaries epe
def tiftosherpagrid(pollutant, variable, sector, aty, net, arr):
    """
     Open Marcos tif and regrid them according to Sherpa's grid
     Input:
     - pollutant, variable, sector, aty, net: to define which file to open
     Output:
     - arr: 2D array with the emissions per sherpa gridd cell

     NB: marco's grid cells are squared, sherpa's grid cell are rectangular
     values in marco's cells are in pairs (in line, 2by2 they have the same
     values so in Sherpa value in a cell is the sum of two adjacent cells
     in Marcos files)

    @author: peduzem
    """
    gdal.UseExceptions()
    ds = None
    try:
        if pollutant != 'CO2':
            ds = gdal.Open('{}_{}/7km_eur_{}_{}{}.tif.tif'.format(pollutant,
                           variable, sector, aty, net))
        else:
            ds = gdal.Open('{}_{}/7km_eur_{}_{}_M{}.tif.tif'.format(pollutant,
                           variable, sector, aty, net))
        em = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
        # re-arrange emission inventory file for Sherpa's grid
        # initialize array
        Amine = em  # Amine : matrix
        for i in range(0, 382):
            ind1 = 2 * (i-1)  # included
            ind2 = 1 + 2 * (i - 1) + 1  # excluded
            Amine[:, i-1] = (np.sum(em[:, ind1:ind2], axis=1))
        Amine[:, 382:384] = 0
        # Cancelling the extra columns and extra rows
        # (is there a better way to do this?**)
        for deli in range(0, 144):
            Amine = np.delete(Amine, (0), axis=0)  # delete first 144 rows
        for delj in range(0, 398):  # delete last 398 columns
            Amine = np.delete(Amine, (383), axis=1)
        Amine_T = Amine[np.newaxis]
        Afinal = np.fliplr(Amine_T)
        arr = Afinal
        return arr
    except(RuntimeError, AttributeError):
        pass
# -----------------------------------------------------------------------------


def readinventories(m_sector, pollutant_lst, sector_lst, net, acttype_lst,
                    nonfuels_lst, name_emi, name_gainsCO2, name_act_fossil,
                    name_act_all, name_ser):
    """
    Read Marco's emission and activity inventory and calculate the reference
    emission factors, for each pollutant, sector, fuel combination

        m_sector = 7 # macro sector, SNAP # 7: road transport
        pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10', 'CO2'] # List of
        pollutants
        sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T',
        'TRA_RD_HDT','TRA_RD_M4' ]# List of subsectors of the m_sector
        net_lst = ['rur', 'urb', 'mot']  ('all')
        acttype_lst = ['GSL', 'MD', 'GAS', 'LPG','TYRE', 'ABRASION','BRAKE']
        # activity types
        nonfuels_lst = ['TYRE', 'ABRASION','BRAKE'] # activities that are not
        fuels
        name_emi= 'emi' name of python binary file to save results
        name_gainsCO2 = 'gainsCO2' name of binary file to save results
        name_ser= 'ser' name of python binary file to save results
        name_act_fossil = 'act_fossil' name of python binary file to save
        results

    OUTPUT:
        act_fossil = activity of fossil fuels [PJ] # taken from the activity
        related to CO2 emissions.
        act_all = activity including fossil and biofuels [PJ] # taken from
        activity realated to PM25 emissions.
        ser = service in [Mvkm]
        emi = emissions of all pollutants [unit]
        gainsCO2 = CO2 emission inventory as in GAINS

    @author: peduzem
    """

    emfco2ipcc = {}
    emfco2ipcc['GSL'] = 69200  # ton/PJ, IPCC2006 69200
    emfco2ipcc['MD'] = 74000  # ton/PJ, IPCC2006 69200
    emfco2ipcc['LPG'] = 63000  # ton/PJ, IPCC2006 69200
    emfco2ipcc['GAS'] = 56100  # ton/PJ, IPCC2006 69200

    # m_sector = 7 # macro sector, SNAP # 7: road transport
    # List of pollutants of interest
    pollutant_lst = ['NH3', 'NOx', 'VOC', 'SO2', 'PM10', 'CO2']
    # List of subsectors of the m_sector
    sector_lst = ['TRA_RD_LD4C', 'TRA_RD_HDB', 'TRA_RD_LD2', 'TRA_RD_LD4T',
                  'TRA_RD_HDT', 'TRA_RD_M4']
    # network (** another option is 'all' but can be extended)
    net_lst = ['rur', 'urb', 'mot'] # 'rur', 'urb', 'mot' or 'all'
#    net_lst =['all']
    # list of activity types
    acttype_lst = ['GSL', 'MD', 'GAS','LPG', 'TYRE', 'ABRASION','BRAKE']
    nonfuels_lst = ['TYRE', 'ABRASION','BRAKE']

    name_emi = 'emi'
    name_gainsCO2 = 'gainsCO2'
    name_ser = 'ser'
    name_act_fossil ='act_fossil'
    name_act_all = 'act_all'

    pollutant = 'CO2'
    variable = 'act_lev'
    act_fossil = {}
    for sector in sector_lst:
        act_fossil[sector]={}
        for aty in [aty for aty in acttype_lst if aty not in nonfuels_lst]:
            act_fossil[sector][aty]={}
            for net in net_lst:
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_fossil)
                if value is not None:
                    act_fossil[sector][aty][net]={}
                    act_fossil[sector][aty][net]=value
            if not act_fossil[sector][aty]:
                act_fossil[sector].pop(aty, None)


    pollutant = 'PM25'
    variable = 'act_lev'
    act_all={}
    for sector in sector_lst:
        act_all[sector]={}
        for aty in [aty for aty in acttype_lst if aty not in nonfuels_lst]:
            act_all[sector][aty]={}
            for net in net_lst:
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_all)
                if value is not None:
                    act_all[sector][aty][net]={}
                    act_all[sector][aty][net]=value
            if not act_all[sector][aty]:
                act_all[sector].pop(aty, None)


    pollutant = 'PM25'
    variable = 'act_lev'
    ser={}
    for sector in sector_lst:
        ser[sector]={}
        for aty in nonfuels_lst:
            ser[sector][aty]={}
            for net in net_lst:
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_all)
                if value is not None:
                    ser[sector][aty][net]={}
                    ser[sector][aty][net]=value
            if not ser[sector][aty]:
                ser[sector].pop(aty, None)

    variable = 'emiss'
    emi = {}
    gainsCO2 = {}
    for pollutant in pollutant_lst:
        emi[pollutant]={}
        for sector in sector_lst:
            emi[pollutant][sector]={}
            gainsCO2[sector]={}
            for aty in acttype_lst:
                emi[pollutant][sector][aty]={}
                gainsCO2[sector][aty]={}
                for net in net_lst:
                    if pollutant=='CO2' and sector in act_fossil and aty in act_fossil[sector]:
                        emi[pollutant][sector][aty][net]={};
                        emi[pollutant][sector][aty][net]=act_fossil[sector][aty][net]*emfco2ipcc[aty]
                        value=tiftosherpagrid(pollutant, variable, sector, aty, net, emi)
                        if value is not None:
                            gainsCO2[sector][aty][net]={}
                            gainsCO2[sector][aty][net]=value
                    elif pollutant != 'CO2':
                        value=tiftosherpagrid(pollutant, variable, sector, aty, net, emi)
                        if value is not None:
                            emi[pollutant][sector][aty][net]={}
                            emi[pollutant][sector][aty][net]=value * 1000 # ton (inventory is in kton)
                if not emi[pollutant][sector][aty]:
                    emi[pollutant][sector][aty].pop(net, None)
                    emi[pollutant][sector].pop(aty, None)

    save_obj((emi), name_emi)
    save_obj((gainsCO2), name_gainsCO2)
    save_obj((act_fossil), name_act_fossil)
    save_obj((act_all), name_act_all)
    save_obj((ser), name_ser)

    return

def aqmeasure(emi, gainsCO2, ser, act_all, act_fossil, area_sel, path_reduction_txt):
    """
    implement aqmeasure TODO
    """
    # calculate the emission factor for every cell in the grid
    # Input parameters
    m_sector = 7 # macro sector, SNAP # 7: road transport
    m_sectorlist = np.arange(1,11)
    # List of pollutants of interest
    pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10', 'CO2']
    # List of subsectors of the m_sector
    sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T','TRA_RD_HDT','TRA_RD_M4' ]
    nonfuels_lst = ['TYRE', 'ABRASION','BRAKE']

    # calculate the implied emission factors, dividing the gridded emissions
    # by the gridded activities (only for the cells where the activity is not
    # zero) TODO: there should be a better way to do this over the stacked dic.
    ef_inv = {}
    for pollutant in emi:
        for sector in emi[pollutant]:
            for aty in emi[pollutant][sector]:
                if aty not in nonfuels_lst:
                    for net in emi[pollutant][sector][aty]:
                        if pollutant=='CO2':
                            ef_inv[pollutant,sector,aty,net] = np.divide(emi[pollutant][sector][aty][net], act_fossil[sector][aty][net], out=np.zeros_like(emi[pollutant][sector][aty][net]), where=act_fossil[sector][aty][net]!=0)
                                  #np.where(act_fossil[sector][aty][net] == 0.0, 0.0, (emi[pollutant][sector][aty][net]/act_fossil[sector][aty][net]))
                        else:
                            ef_inv[pollutant,sector,aty,net] = np.divide(emi[pollutant][sector][aty][net], act_all[sector][aty][net], out=np.zeros_like(emi[pollutant][sector][aty][net]), where=act_all[sector][aty][net]!=0)
                            #= np.where(act_all[sector][aty][net] == 0.0, 0.0, (emi[pollutant][sector][aty][net]/act_all[sector][aty][net]))
                if aty in nonfuels_lst:
                    for net in emi[pollutant][sector][aty]:
                        ef_inv[pollutant,sector,aty,net] = np.divide(emi[pollutant][sector][aty][net], ser[sector]['TYRE'][net], out=np.zeros_like(emi[pollutant][sector][aty][net]), where=ser[sector]['TYRE'][net]!=0)
                              #np.where(ser[sector]['TYRE'][net] == 0.0, 0.0, (emi[pollutant][sector][aty][net]/ser[sector]['TYRE'][net]))

    # technical measures
    # New emission factors after the implementation of measures:
    # first option - reduction of the emission factors for each sector/activity
    # read csv file with the emission factors
    df = pd.read_csv('input/ef_red_sherpaeco.csv',  index_col=[0,1], names=['POLL','ACT'].extend(sector_lst), skipinitialspace=True)
    # second option (to be implemented**)- emission factors for the best available technology
    # read the precursors list (as in model1)

    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    # create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)

   # caluclate total emission in the area of interest for the base case (Sherpa's emission inventory)
    em_bc ={} # emissions for the base case
    em_bc_t_percell = {}
    for precursor in precursor_lst:
        for m_sec in m_sectorlist:
            em_bc[precursor,m_sec-1] = np.sum(emission_dict[precursor][m_sec-1] * area * (area_sel/100))
        em_bc_t_percell[precursor]= np.sum(emission_dict[precursor][m_sec-1] * area * area_sel/100 for m_sec in m_sectorlist)

    # calculate delta emissions and inventory emissions per cell (in the area of interest)
    delta_em_percell={}
    em_percell={}
    for pollutant in emi:
        for sector in emi[pollutant]:
           for aty in emi[pollutant][sector]:
                for net in emi[pollutant][sector][aty]:
                   # ef = ef_inv[pollutant,sector,aty,net] *(1 - df.loc[pollutant,aty][sector])
                    if aty not in nonfuels_lst:
                        if pollutant != 'CO2':
                            delta_em_percell[pollutant,sector,aty,net]=act_all[sector][aty][net]*ef_inv[pollutant,sector,aty,net]*df.loc[pollutant,aty][sector]* (area_sel/100)
                            em_percell[pollutant,sector,aty,net]=act_all[sector][aty][net]*ef_inv[pollutant,sector,aty,net] * (area_sel/100)
                        if pollutant=='CO2':
                            delta_em_percell[pollutant,sector,aty,net]=act_fossil[sector][aty][net]*ef_inv[pollutant,sector,aty,net]*df.loc[pollutant,aty][sector]*area_sel/100
                            em_percell[pollutant,sector,aty,net]=act_fossil[sector][aty][net]*ef_inv[pollutant,sector,aty,net] * (area_sel/100)
                    if aty in nonfuels_lst:
                            delta_em_percell[pollutant,sector,aty,net]=ser[sector][aty][net]*ef_inv[pollutant,sector,aty,net]*df.loc[pollutant,aty][sector] * (area_sel/100)
                            em_percell[pollutant,sector,aty,net]=ser[sector][aty][net]*ef_inv[pollutant,sector,aty,net] * (area_sel/100)


    # calculate delta emissions per pollutant and emissions per pollutant per cell and total from Marcos inventory
    delta_em_pp_percell={}
    em_pp_percell={}
    delta_em_pp_t={}
    em_pp_t={}
    for pollutant in pollutant_lst:
        delta_em_pp_percell[pollutant] = np.sum(np.sum(np.sum(delta_em_percell[pollutant,sector,aty,net] for net in emi[pollutant][sector][aty]) for aty in emi[pollutant][sector]) for sector in emi[pollutant])
        em_pp_percell[pollutant] = np.sum(np.sum(np.sum(em_percell[pollutant,sector,aty,net] for net in emi[pollutant][sector][aty]) for aty in emi[pollutant][sector]) for sector in emi[pollutant])
        delta_em_pp_t[pollutant]=np.sum(delta_em_pp_percell[pollutant])
        em_pp_t[pollutant]=np.sum(em_pp_percell[pollutant])
        print(pollutant, delta_em_pp_t[pollutant])

    delta_emission = {}
    emissions = {}
    emissions_percell = {}
    for precursor in precursor_lst:
        for pollutant in pollutant_lst:
            if pollutant == precursor:
                delta_emission[precursor] = delta_em_pp_percell[pollutant]
                emissions[precursor]=em_pp_t[pollutant]
                emissions_percell[precursor]=em_pp_percell[pollutant]
            if precursor == 'NMVOC':
                delta_emission[precursor] = (delta_em_pp_percell['VOC'])
                emissions[precursor]=em_pp_t['VOC']
                emissions_percell[precursor]=em_pp_percell['VOC']
            if precursor == 'SOx':
                delta_emission[precursor] = (delta_em_pp_percell['SO2'])
                emissions[precursor]=em_pp_t['SO2']
                emissions_percell[precursor]=em_pp_percell['SO2']
            if precursor == 'PPM':
                delta_emission[precursor] = (delta_em_pp_percell['PM10'])
                emissions[precursor]=em_pp_t['PM10']
                emissions_percell[precursor]=em_pp_percell['PM10']

#   Reductions per precursor
    red = {}
    m_sector = 7
    for precursor in precursor_lst:
        red[precursor]={}
        red[precursor] = np.sum(delta_emission[precursor])/em_bc[precursor,m_sector-1]*100
#        red[precursor] = np.sum(delta_emission[precursor])/emissions[precursor]*100

    delta_emission_dict={}
    for precursor in precursor_lst:
#        delta_emission_dict[precursor] = np.sum((delta_emission[precursor]/area), axis=0) #** check this (this is to sum over the macrosector whe I will have that information)
        delta_emission_dict[precursor] = np.sum((delta_emission[precursor]/area), axis=0) #delta_emission[precursor]/area

#    # constant over the whole area... I think this is what GAINS does. TODO has to be in percentage not absolute value!!
    delta_emission_dict_per={}
    for precursor in precursor_lst:
#        delta_emission_dict_per[precursor] = (np.sum(delta_emission[precursor])*area_sel/(np.sum(area_sel[(0)])))/area
         delta_emission_dict_per[precursor] = (em_bc_t_percell[precursor] * np.sum(delta_emission[precursor]) / np.sum(em_bc_t_percell[precursor]))/area
#                                (np.sum(delta_emission[precursor])*area_sel[(0)]/(np.sum(area_sel[(0)])))/area

    reductions = {}
    m_sector = 7
    reductions[m_sector-1]={}
#                                   NOx	        NMVOC       NH3           PPM        SOx
    reductions[m_sector-1] = [red['NOx'], red['NMVOC'], red['NH3'] , red['PPM'], red['SOx']]

    write_reductions(path_reduction_txt, reductions[m_sector-1])

    return delta_emission_dict, delta_emission_dict_per

# -------------------------------------------------------------------------
# main program starts here
# -------------------------------------------------------------------------

if __name__ == '__main__':
    import os, errno
    from postcompute_healthia import main_healthimpact
    # -------------------------------------------------------------------------
    # Preparing the input data
    # -------------------------------------------------------------------------
    workdir = 'workdir\\sherpaeco\\'
    # Open model file, I just need it to load lons and lats
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(lon_array)
    n_lat = len(lat_array)
    rootgrp.close()

    # Area of each cell in the domain
    nc_area = 'input/JRC01.nc'
    rootgrp = Dataset(nc_area, mode='r')
    area = rootgrp.variables['surface'][:]
    area_units = rootgrp.variables['surface'].units
    rootgrp.close()

    # info on areas and percentage of grids in areas
    grid_txt = 'input\\selection\\grid_intersect'
    gcities_txt = 'input\\selection\\grid_int_gcities'
    fua_txt = 'input\\selection\\\\grid_int_fua'


    path_reduction_txt='input/sherpaeco_reduction.txt'




    # -------------------------------------------------------------------------

    # uncomment to read inventories from Marco
    m_sector = 7 # macro sector, SNAP # 7: road transport
    m_sectorlist = np.arange(1,11)
    # List of pollutants of interest
#    pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10', 'CO2']
#    # List of subsectors of the m_sector
#    sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T','TRA_RD_HDT','TRA_RD_M4' ]
#    nonfuels_lst = ['TYRE', 'ABRASION','BRAKE']
#    readinventories(m_sector, pollutant_lst, sector_lst, net_lst, acttype_lst,
#                    nonfuels_lst, name_emi, name_gainsCO2, name_act_fossil,
#                    name_act_all,name_ser)

    # -------------------------------------------------------------------------
    # Measures implementation (**will be translated to a def)
    # -------------------------------------------------------------------------
    eu28 = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES','FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT','NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
#    eu28 = ['DE']
    # uncomment to implement measure TODO
    name_emi = 'emi'
    name_gainsCO2 = 'gainsCO2'
    name_ser = 'ser'
    name_act_fossil = 'act_fossil'
    name_act_all = 'act_all'
    # load results of readinventories
    emi = load_obj(workdir + name_emi)
    gainsCO2 = load_obj(workdir + name_gainsCO2)
    ser = load_obj(workdir + name_ser)
    act_all = load_obj(workdir + name_act_all)
    act_fossil = load_obj(workdir + name_act_fossil)

    # Target area:
    level = 'NUTS_Lv0'
    parea = 'parea'
    
    deltayll_reg = {}
    delta_mort_reg = {}
    deltayll_spec_reg = {}
    nc_selarea = 'workdir\\area_IT.nc'
    rootgrp = Dataset(nc_selarea, mode='r')
    area_sel = rootgrp.variables['AREA'][:]
    # area_units = rootgrp.variables['AREA'].units
    rootgrp.close()
    eu28_area = np.zeros(np.shape(area_sel))
#    for code in eu28:
#
#        # path to store the selected area
#    #    nc_selarea = 'workdir\\area_{}{}.nc'.format(level, code)
#        nc_selarea = 'workdir\\area_{}.nc'.format(code)        
#        # reduction area
#        
#        if not os.path.exists(nc_selarea):
#            area_sel = gridint_toarray(level, 'parea', code)
#            write_nc(area_sel, nc_selarea, 'AREA', '%')
#        else: 
#            rootgrp = Dataset(nc_selarea, mode='r')
#            area_sel = rootgrp.variables['AREA'][:]
#            # area_units = rootgrp.variables['AREA'].units
#            rootgrp.close()
#            
#        eu28_area = eu28_area + area_sel
#    
#    nc_selareaeu28 = 'workdir\\area_EU28.nc'   
#    if not os.path.exists(nc_selarea):
#            area_sel = gridint_toarray(level, 'parea', code)
#            write_nc(eu28_area, nc_selareaeu28, 'AREA', '%')
##    write_nc(eu28_area, nc_selareaeu28, 'AREA', '%')
#    
#    delta_emission_dict, delta_emission_dict_per = aqmeasure(
#                emi, gainsCO2, ser, act_all, act_fossil, eu28_area,
#                path_reduction_txt)
#    
#        # -------------------------------------------------------------------------
#        # Running module1
#        # -------------------------------------------------------------------------
#        # if it doesn't exist start=0 and dividsor=1
#    progresslog = 'input/progress.log'
#    # run module 1 with progress log
#    output = 'outputdieselgate//'
#
#    try:
#        os.makedirs(output)
#    except OSError as e:
#        if e.errno != errno.EEXIST:
#            raise
#
#    proglog_filename = path_result_cdf_test + 'proglog'
#    write_progress_log(proglog_filename, 25, 2)
#    start = time()
#    shrp.module1(path_emission_cdf_test, nc_selareaeu28,
#                 path_reduction_txt, path_base_conc_cdf_test,
#                 path_model_cdf_test, output)
#
#    stop = time()
#    print('Module 1 run time: %s sec.' % (stop-start))
#    remove(proglog_filename)
#    
    
    output = 'outputdieselgate//percountry//'

    try:
        os.makedirs(output)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    countries =['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES','FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT','NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK'] 
#    countries = ['DE']
    for code in countries:
        outputpc = output + code + '//'
        try:
            os.makedirs(outputpc)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # path to store the selected area
    #    nc_selarea = 'workdir\\area_{}{}.nc'.format(level, code)
        nc_selarea = 'workdir//area_{}.nc'.format(code)        
        # reduction area
        
        if not os.path.exists(nc_selarea):
            area_sel = gridint_toarray(level, 'parea', code)
            write_nc(area_sel, nc_selarea, 'AREA', '%')
        else: 
            rootgrp = Dataset(nc_selarea, mode='r')
            area_sel = rootgrp.variables['AREA'][:]
            # area_units = rootgrp.variables['AREA'].units
            rootgrp.close()
            delta_emission_dict, delta_emission_dict_per = aqmeasure(
                emi, gainsCO2, ser, act_all, act_fossil, area_sel,
                path_reduction_txt)
        proglog_filename = path_result_cdf_test + 'proglog'
        write_progress_log(proglog_filename, 25, 2)
        start = time()
        shrp.module1(path_emission_cdf_test, nc_selarea,
                 path_reduction_txt, path_base_conc_cdf_test,
                 path_model_cdf_test, outputpc)

        stop = time()
        print('Module 1 run time: %s sec.' % (stop-start))
        remove(proglog_filename)
    
    #    progresslog = 'input/progress.log'
    #    # run module 1 with progress log
    #    output = 'output2/'
    #    proglog_filename = path_result_cdf_test + 'proglog'
    #    write_progress_log(proglog_filename, 25, 2)
    #    start = time()
    #    shrp1.module1(path_emission_cdf_test, nc_selarea,
    #                  delta_emission_dict, path_base_conc_cdf_test,
    #                  path_model_cdf_test, output)
    #    stop = time()
    #    print('Module 1 run time: %s sec.' % (stop-start))
    #    remove(proglog_filename)
    #
    #    output = 'output3/'
    #    proglog_filename = path_result_cdf_test + 'proglog'
    #    write_progress_log(proglog_filename, 25, 2)
    #    start = time()
    #    shrp1.module1(path_emission_cdf_test, nc_selarea,
    #                  delta_emission_dict_per, path_base_conc_cdf_test,
    #                  path_model_cdf_test, output)
    #    stop = time()
    #    print('Module 1 run time: %s sec.' % (stop-start))
    #    remove(proglog_filename)
    
    
         # -------------------------------------------------------------------------
         # (Cost) Benefit Analysis
         # -------------------------------------------------------------------------
    
    deltaconc = 'outputdieselgate/delta_concentration.nc'

#        deltayll_reg[code], delta_mort_reg[code], deltayll_spec_reg[code] = calc_impacts(
#                deltaconc, nc_selarea, path_tiff, path_mortbaseline,
#                path_dust_conc_cdf_test, path_salt_conc_cdf_test,
#                code)

#    deltaconc2 = 'output2/delta_concentration.nc'
#    deltayll_reg2, delta_mort_reg2, deltayll_spec_reg2 = calc_impacts(
#            deltaconc2, nc_selarea, path_tiff, path_mortbaseline,
#            path_dust_conc_cdf_test, path_salt_conc_cdf_test,
#            code)
#    
#    deltaconc3 = 'output3/delta_concentration.nc'
#    deltayll_reg3, delta_mort_reg3, deltayll_spec_reg3 = calc_impacts(
#            deltaconc3, nc_selarea, path_tiff, path_mortbaseline,
#            path_dust_conc_cdf_test, path_salt_conc_cdf_test,
#            code)    
    
    
#     -------------------------------------------------------------------------
    # Output of results (todo)
    # -------------------------------------------------------------------------

##    plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
#    # create new netcdf file for results
#    nc_out = 'output/delta_yll.nc'
#    fh = Dataset(nc_out, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude',  n_lat)
#    fh.createDimension('longitude', n_lon) #x
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    yllout = fh.createVariable('yll', 'f4', ('latitude', 'longitude'))
#    fh.variables['yll'].units = 'years'
#    fh.variables['yll'].long_name = 'years of life loss'
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    yllout[:] = delta_yll3 #area_sel[(0)]*100
#    fh.close()

#    # create new netcdf file for results
#    nc_out = 'output2/delta_yll.nc'
#    fh = Dataset(nc_out, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude',  n_lat)
#    fh.createDimension('longitude', n_lon) #x
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    yllout = fh.createVariable('yll', 'f4', ('latitude', 'longitude'))
#    fh.variables['yll'].units = 'years'
#    fh.variables['yll'].long_name = 'years of life loss'
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    yllout[:] = delta_yll2
#    fh.close()
#
#        # create new netcdf file for results
#    nc_out = 'output3/delta_yll.nc'
#    fh = Dataset(nc_out, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude',  n_lat)
#    fh.createDimension('longitude', n_lon) #x
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    yllout = fh.createVariable('yll', 'f4', ('latitude', 'longitude'))
#    fh.variables['yll'].units = 'years'
#    fh.variables['yll'].long_name = 'years of life loss'
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    yllout[:] = delta_yll3
#    fh.close()
#

# -------------------------------------------------------------------------
#    # Check inventory
#    # initialization
#    driver = ogr.GetDriverByName('ESRI Shapefile')
#    shp = driver.Open(r'\shapefiles')
#    layer = shp.GetLayer()
#
#    rootgrp = Dataset(path_model_cdf_test, 'r')
#    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
#    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)
#    em_dict_area_t = {}
#    for precursor in precursor_lst:
#        for snap in range(1, 11):
#            em_dict_area_t[precursor, snap - 1] = np.sum(
#            emission_dict[precursor][snap - 1]  * area * area_sel[(0)])  # [Mg]
#
#    ser_area_t = {} # total service [Mvkm] in the area of interest
#    act_all_area_t = {} # total activity [PJ] in the area of interest
#    act_foss_area_t = {} # total activity [PJ] in the area of interest
#    emiCO2_area_t = {}
#    emiPM10_area_t={}
#    emiNOx_area_t={}
#    gainsCO2_area_t = {}
#    emiNH3_area_t = {}
#    emiSO2_area_t ={}
#
#
#    for sector in act_all:
#        for aty in act_all[sector]:
#            for net in act_all[sector][aty]:
##                act_all_area[sector,aty,net] = act_all[sector][aty][net] * area_area / 100  # [PJ]
#                act_all_area_t[sector,aty,net] = np.sum(act_all[sector][aty][net] * area_sel[(0)] )  # [PJ]
#                act_foss_area_t[sector,aty,net] = np.sum(act_fossil[sector][aty][net] * area_sel[(0)])  # [PJ]
#
##                for pollutant in ['CO2']:
##                    try:
#                emiCO2_area_t[sector, aty, net]=np.sum(emi['CO2'][sector][aty][net] * area_sel[(0)])
#                gainsCO2_area_t[sector, aty, net]=np.sum(gainsCO2[sector][aty][net] * area_sel[(0)])
#                emiPM10_area_t[sector, aty, net]=np.sum(emi['PM10'][sector][aty][net] * area_sel[(0)])
#                emiNOx_area_t[sector, aty, net]=np.sum(emi['NOx'][sector][aty][net] * area_sel[(0)])
#                emiNH3_area_t[sector, aty, net]=np.sum(emi['NH3'][sector][aty][net] * area_sel[(0)])
#                try:
#                    emiSO2_area_t[sector, aty, net]=np.sum(emi['SO2'][sector][aty][net] * area_sel[(0)])
#                except(KeyError):
#                    emiSO2_area_t[sector, aty, net]=0
#    pd.set_option('precision', 5)
##   act_fossil_area_t[sector,aty,net] = np.sum(act_fossil_area[sector,aty,net])  # [PJ]
#    # create dataframe with activities (all, fossil) and emissions
#    data= np.transpose([(list(act_all_area_t.values())),
#                        (list(act_foss_area_t.values())),
##                        (list(gainsCO2_area_t.values())),
##                        (list(emiPM10_area_t.values())),
##                        (list(emiNOx_area_t.values())),
##                        (list(emiSO2_area_t.values())),
#                        (list(emiNH3_area_t.values()))])
#    index = pd.MultiIndex.from_tuples(list(act_all_area_t.keys()), names=['sector', 'aty','net'])
##    dfact=pd.DataFrame(data,index, columns =['act_all', 'act_foss','CO2','PM10','NOx','SO2', 'NH3']).sortlevel()
#    dfact=pd.DataFrame(data,index, columns =['act_all', 'act_foss','NH3']).sortlevel()
#
#    comparison=dfact.groupby(level=1).sum()
#    grouped=dfact.groupby(level=['sector', 'aty']).sum()
#    grouped2=dfact.groupby(level=['sector', 'aty']).apply(lambda x:
#                                                 100 * x / x.sum())


#    state_pcts = state_office.groupby(level=0).apply(lambda x:
#                                                 100 * x / float(x.sum()))

# -------------------------------------------------------------------------

#    act_all_area_t_net={}
#    for sector in act_all:
#            for aty in act_all[sector]:
#                act_all_area_t_net[sector,aty]=np.sum(act_all_area_t[sector,aty,net] for net in act_all[sector][aty])


#    act_all_pp_area_t={}
#    for net in emi['CO2'][sector][aty]:
#        act_all_pp_area_t[net] = np.sum(np.sum(act_all_area_t[sector,aty,net]  for aty in emi['CO2'][sector]) for sector in emi['CO2'])
#
#    emi_area_t ={}
#    emi_area ={}
#    for pollutant in emi:
#        for sector in emi[pollutant]:
#            for aty in emi[pollutant][sector]:
#                for net in emi[pollutant][sector][aty]:
#                    emi_area[pollutant,sector,aty,net] = emi[pollutant][sector][aty][net] * area_area / 100
#                    emi_area_t[pollutant,sector,aty,net] = np.sum(emi_area[pollutant,sector,aty,net])

#    em_all_pp_area_t_net={}
#    for net in emi['CO2'][sector][aty]:
#        em_all_pp_area_t_net[net] = np.sum(np.sum(emi_area_t['VOC',sector,aty,net]  for aty in emi['CO2'][sector]) for sector in emi['CO2'])

#    ser_area_t={}
#    ser_area={}
#    ser_all = {}
#    for sector in ser:
#        for net in ser[sector]['TYRE']:
#            ser_all[sector,net] = ser[sector]['TYRE'][net]
#            ser_area[sector,net] = ser[sector]['TYRE'][net] * area_area / 100 # [Mvkm]
#            ser_area_t[sector,net] = np.sum(ser_area[sector,net]) # [Mvkm]

#    ef_area_t = {}
#    ef_area = {}
##    for sec in sector_lst: dct[sec] = {}
#    [(x, y) for x in [1,2,3] for y in [3,1,4] if x != y]
##    dict([('sape', 4139), ('guido', 4127), ('jack', 4098)])
##    dct = defaultdict(dict)
#
#    dct1 = {}
## dict(pol,{})
#     #    dct[pollutant] = {} for pollutant in pollutant_lst
#    dct1 =[dict([((pol,sector,aty,net), {})])
#        for pol in pollutant_lst for sector in sector_lst
#        for aty in acttype_lst for net in net_lst]
#
#    dct2 = {}
##    for pol in pollutant_lst:
##        for sector in sector_lst:
##            for aty in acttype_lst:
##                for net in net_lst:
##                    dct2[pol,sector,aty,net] ={}
##setdefault(
##    d = dict()
##    d =dict(d.setdefault(pol, {}).setdefault(sector, {}).setdefault(aty, {}) for pol in pollutant_lst
##        for sector in sector_lst
##        for aty in acttype_lst)
#
##     dct.update(dict(pol,{}) for pol in pollutant_lst)
##    dct = [dict(pol, {}) for pol in pollutant_lst]
##    dct['CO2']=3
##    em ={} # emission after measure implementation
##    for pollutant in pollutant_lst:
##        em[pollutant]={}
#
##    target_dict = defaultdict(dict)
##    target_dict['PM10']['TRA_RD_HDT'] = 5
#
##    [aty for aty in acttype_lst if aty not in nonfuels_lst]


    # -------------------------------------------------------------------------
    # check grid (marcos vs sherpas)
    # -------------------------------------------------------------------------

#    rootgrp = Dataset('7km_eur_TRA_RD_LD4C_GSL_Mall.nc', 'r')
#    lon_marconc = rootgrp.variables['lon'][:]
#    lat_marconc = rootgrp.variables['lat'][:]
#    ncem=rootgrp.variables['7km_eur_TRA_RD_LD4C_GSL_Mall.tif.tif'][:]
#    rootgrp.close()
#
#    ds   = gdal.Open('CO2_emiss/7km_eur_TRA_RD_LD4C_GSL_Mall.tif.tif')
#    arr    = ds.ReadAsArray()
#    [cols,rows] = arr.shape
#    (Xarr, deltaX, rotation, Yarr, rotation, deltaY) = ds.GetGeoTransform()
#    CO2_TRA_RD_LD4C_GSL_Mall = np.array(ds.GetRasterBand(1).ReadAsArray())
#    ds = None
#    emidct={}
#    # Define longitude and latitude for the data from marco
#    longl = []
#    for i in range(0, rows):
#        longl.append(Xarr + i*deltaX + deltaX*0.5)
#    lonsmtif=np.array(longl)
#    latl = []
#    for i in range(0, cols):
#        latl.append(Yarr + i*deltaY + deltaY*0.5)
#    latsmtif=np.array(latl)
#    X, Y = np.meshgrid(lonsmtif , latsmtif)
#
#    Amine = CO2_TRA_RD_LD4C_GSL_Mall # Amine : matrix
#    for i in range(0,382):
#        ind1=2*(i-1)  # included
#        ind2=1+2*(i-1)+1 # excluded
#        Amine[:,i-1]=(np.sum(CO2_TRA_RD_LD4C_GSL_Mall[:,ind1:ind2],axis=1))
#    Amine[:,382:384]=0
#    # Cancelling the extra columns and extra rows
#    # (there has to be a better way to do this)
#    for deli in range (0,144):
#        Amine = np.delete(Amine, (0), axis=0) # delete first 144 rows
#    for delj in range (0,398): # delete last 398 columns
#        Amine = np.delete(Amine, (383), axis=1)
#    Amine_T = Amine[np.newaxis]
#    Afinal=np.fliplr(Amine_T)
#    emidct=Afinal
#
#
#
    # create new netcdf file for results
#    nc_file3 = 'netcdf/myTRA_RD_LD4C_GSL_Mall.nc'
#    co2em = Dataset(nc_file3, mode='w', format='NETCDF3_CLASSIC')
#    co2em.createDimension('latitude',  len(emission_dict['lat_array']))
#    co2em.createDimension('longitude', len(emission_dict['lon_array'])) #x
#    latitude = co2em.createVariable('latitude', 'f4', ('latitude',))
#    longitude = co2em.createVariable('longitude', 'f4', ('longitude',))
#    emissions = co2em.createVariable('emissions', 'f4', ('latitude', 'longitude'))
#    co2em.variables['emissions'].units = 'unit'
#    co2em.variables['emissions'].long_name = 'unit'
#    longitude[:] = emission_dict['lon_array']
#    latitude[:] = emission_dict['lat_array']
#    emissions[:] = ef_area['CO2','TRA_RD_LD4T','GSL','all']
#    co2em.close()
##
#    nc_file3 = 'netcdf/trial.nc'
#    co2em = Dataset(nc_file3, mode='w', format='NETCDF3_CLASSIC')
#    co2em.createDimension('latitude',  len(emission_dict['lat_array']))
#    co2em.createDimension('longitude', len(emission_dict['lon_array'])) #x
#    latitude = co2em.createVariable('latitude', 'f4', ('latitude',))
#    longitude = co2em.createVariable('longitude', 'f4', ('longitude',))
#    emissions = co2em.createVariable('emissions', 'f4', ('latitude', 'longitude'))
#    co2em.variables['emissions'].units = 'unit'
#    co2em.variables['emissions'].long_name = 'unit'
#    longitude[:] = emission_dict['lon_array']
#    latitude[:] = emission_dict['lat_array']
#    emissions[:] = area_area
##    emissions[:] = np.where(act_all['TRA_RD_LD4C']['GSL']['mot']==0, 0, (((emi['VOC']['TRA_RD_LD4C']['GSL']['mot'])/(act_all['TRA_RD_4C']['GSL']['mot']))))
#    co2em.close()
#
##
#    nc_file3 = 'netcdf/deltayllpploriginal.nc'
#    co2em = Dataset(nc_file3, mode='w', format='NETCDF3_CLASSIC')
#    co2em.createDimension('latitude',  len(emission_dict['lat_array']))
#    co2em.createDimension('longitude', len(emission_dict['lon_array'])) #x
#    latitude = co2em.createVariable('latitude', 'f4', ('latitude',))
#    longitude = co2em.createVariable('longitude', 'f4', ('longitude',))
#    emissions = co2em.createVariable('emissions', 'f4', ('latitude', 'longitude'))
#    co2em.variables['emissions'].units = 'unit'
#    co2em.variables['emissions'].long_name = 'unit'
#    longitude[:] = emission_dict['lon_array']
#    latitude[:] = emission_dict['lat_array']
##    emissions[:] = ef_all['CO2','TRA_RD_LD4T','GSL','mot']
#    emissions[:] = ef_inv['NOx','TRA_RD_HDB','MD','mot'] #np.where(act_fossil['TRA_RD_LD4T']['GSL']['all']==0, 0, (emi['CO2']['TRA_RD_LD4T']['GSL']['all']/act_fossil['TRA_RD_LD4T']['GSL']['all']))
#    co2em.close()

#    np.sum(emi_area['NOx','TRA_RD_LD4C','GSL','urb']-em_new['NOx','TRA_RD_LD4C','GSL','urb'])


#    H = area_sel[(0)]
##
#    fig = plt.figure(figsize=(6, 3.2))
#
#    ax = fig.add_subplot(111)
#    ax.set_title('colorMap')
#    plt.imshow(H)
#    ax.set_aspect('equal')
#
#    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
#    cax.get_xaxis().set_visible(False)
#    cax.get_yaxis().set_visible(False)
#    cax.patch.set_alpha(0)
#    cax.set_frame_on(False)
#    plt.colorbar(orientation='vertical')
#    plt.show()
