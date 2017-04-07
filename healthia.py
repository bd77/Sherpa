# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017

@author: peduzem
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset  # for using netCDF files
import numpy as np
import scipy.io as sio

from module7_custom_ep import read_nuts_area
from sherpa_globals import (path_area_cdf_test, path_model_cdf_test,
                            path_emission_cdf_test, path_result_cdf_test)
from sherpa_auxiliaries import (create_emission_reduction_dict,
                                create_emission_dict, create_window,
                                read_progress_log, write_progress_log)
from sherpaeco import gridint_toarray

def calc_impacts(concncdf, nc_area, spec='d'):

    """
    Health impacts of PM2.5
    Methodology detailed in the REVIHAAP and HRAPIE studies led by WHO-Europe,
    as described in:
    Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios for the
    EU Clean Air Package Version. Version 2

    Concentration-Response-Function taken from the software from AirQ+ (WHO)

    Years of life loss as calculated as a function of life expectancy according to
    Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios for the EU Clean Air Package Version. Version 2

    Data for baseline  population
    ICD codes: ICD-10: A00-B99,C00-D48,D50-D89,E00-E88,F01-F99,G00-G98,
    H00-H59,H60-H93,I00-I99,J00-J98,K00-K92,L00-L98,M00-M99,N00-N98,
    O00-O99,P00-P96,Q00-Q99,R00-R99
    Age: '30 - 85 +'
    Sex: Both
    http://data.euro.who.int/dmdb/
    [Accessed December 13, 2016].

    NB: Hereafter a positive delta means a reduction!


    Input:
        - nc_area: path to netcdf file of area of interest (nc file with 100 in the cells of the area of interest)
        - deltaconc: path to the file  with the delta concentrations (output of module1)
    Output:
        - deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg, delta_yl

    @author: peduzem
    """
    if spec == 'a':
        rootgrp = Dataset(concncdf, mode='r')
        pm25_conc = rootgrp.variables['conc'][:]
        rootgrp.close()
    else:
        fh_deltapm25 = Dataset(concncdf, mode='r')
        pm25_conc = fh_deltapm25.variables['delta_concentration'][:]
#    pm25_delta = fh_deltapm25.variables['conc'][:]
        fh_deltapm25.close()

    # Load population file from Marco (treated in "regridPop" routine from Enrico).
    # The matlab file contains a structure called Anew,
    # I carried out the same operations as in Enrico's matlab file.
    A = sio.loadmat('input/population.mat')
    Anew = A['Anew']
    Anew_T = Anew[np.newaxis]
    popall=np.fliplr(Anew_T)

#    nc_pop='input/population.nc'
#    fh = Dataset(nc_pop, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude',  len(lat_array))
#    fh.createDimension('longitude', len(lon_array)) #x
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    pop = fh.createVariable('pop', 'f8', ('latitude', 'longitude'))
#    fh.variables['pop'].units = 'ppl'
#    fh.variables['pop'].long_name = 'population'
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    pop[:] = popall
#    fh.close()

    rootgrp = Dataset(nc_area, mode='r')
    area_area = rootgrp.variables['AREA'][:]
    rootgrp.close()

    # total pop according to Marco (all age groups)
    sumpop=np.sum(popall)

    # TOTAL population above 30 in EU28 data.euro.who.int/dmdb/
    pop = 331923577
    # TOTAL DEATHS  above 30 in EU28 data.euro.who.int/dmdb/
    deaths = 4639244

#    # Potential Years of life loss, PYLL 30+ total in EU28 data.euro.who.int/dmdb/
#    ylltot = 14038453.71

    # Distribution of the population over 30
    # HP: assuming it has the same distribution as all the population (as provided by Marco)
    pop30plus = (popall/sumpop) * pop

    # open delta concentration for PM25 - result of module1


    # Death rate over 30
    drate = deaths/pop

    # crf derived by RR in Anenberg, S.C. et al., 2010.
    # crf = 0.006

    # crf Taken from AirQ+ (WHO)
#    cutoff = 10 # microg/m3 # Taken from AirQ+ (WHO)
    mbeta = 0.006015392281974714 # Taken from AirQ+ (WHO)
    lbeta = 0.003922071315328133 # Taken from AirQ+ (WHO)
    hbeta = 0.007973496801885352 # Taken from AirQ+ (WHO)
    beta = [lbeta, mbeta, hbeta]
#    baseconfile = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
#    fh_basecon = Dataset(baseconfile , mode='r')
#    baseconc = fh_basecon.variables['conc'][:]
#    fh_basecon.close()
    #pddeltaconc.loc['delta_concentration'].loc['delta_concentration']

    # Delta mortality: uncomment the most suitable one
    # linear approximation
    # delta_mort = pm25_delta*crf*pop30plus*drate
    # exponential approximation
    delta_mort = [(1-(np.exp(-beta[i]*pm25_conc)))*pop30plus*drate for i in range(len(beta))]
    delta_mort_tot = [np.sum(delta_mort[i]) for i in range(len(beta))]
    delta_mort_reg = [np.sum(delta_mort[i]*area_area/100) for i in range(len(beta))]

    # Delta Years of life loss (yll) according to Anenberg 2010
    # delta_yll = delta_mort * ylltot / deaths
    # ** I tried validating these values but I obtain very different results
    # from the ones in Holland, M., 2014.
    # 79.9 is the life expectancy, should be by country, this is the average from EUROSTAT
    # http://ec.europa.eu/eurostat/statistics-explained/images/9/93/Life_expectancy_at_birth%2C_1980%E2%80%932014_%28years%29_YB16.png
    lyg = np.exp(8.161-(0.04478*79.9)) # life years gained per 100000 ppl for a unit concentration change
    delta_yll= pm25_conc * lyg / 100000* pop30plus
    # Calculate delta yll in the selected region
    deltayll_reg = np.sum(delta_yll * area_area/100)
    # yll gained per 100000 people
    deltayll_spec_reg = (
            (deltayll_reg / np.sum(pop30plus * area_area/100)) * 100000)

    # Calculate delta yll in total
    deltayll_tot = np.sum(delta_yll)

#    # create new netcdf file for results
#    nc_out = 'netcdf/pop.nc'
#    fh = Dataset(nc_out, mode='w', format='NETCDF3_CLASSIC')
#    fh.createDimension('latitude',  n_lat)
#    fh.createDimension('longitude', n_lon) #x
#    latitude = fh.createVariable('latitude', 'f4', ('latitude',))
#    longitude = fh.createVariable('longitude', 'f4', ('longitude',))
#    yllout = fh.createVariable('pop', 'f4', ('latitude', 'longitude'))
#    fh.variables['pop'].units = 'people'
#    fh.variables['pop'].long_name = 'number of people'
#    longitude[:] = lon_array
#    latitude[:] = lat_array
#    yllout[:] = pop30plus
#    fh.close()
#
    return (deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg,
            deltayll_spec_reg)

# -------------------------------------------------------------------------
# main program starts here
# -------------------------------------------------------------------------

if __name__ == '__main__':

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

    # Target area:
    level = 'FUA_CODE'
    parea = 'parea'
    code = 'IT002L2'

    # path to store the selected area
    nc_selarea = 'workdir/selarea.nc'

    path_reduction_txt='input/sherpaeco_reduction.txt'

    #reduction area
    area_sel = gridint_toarray(
            level, parea, code, lat_array, lon_array, grid_txt, gcities_txt,
            fua_txt, nc_selarea)

    #deltaconc='output/delta_concentration.nc'
    deltaconc='input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
    deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg, deltayll_spec_reg = calc_impacts(deltaconc, nc_selarea, spec='a')
#
    H=area_sel[(0)]*100
#
    fig = plt.figure(figsize=(6, 3.2))

    ax = fig.add_subplot(111)
    ax.set_title('colorMap')
    plt.imshow(H)
    ax.set_aspect('equal')

    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.show()


