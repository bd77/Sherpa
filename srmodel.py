'''
Created on 8 Feb 2017

@author: degraba
'''

from numpy import zeros, power, sqrt, ones, savetxt, sum
from time import time
from math import isnan, floor
from netCDF4 import Dataset


# auxiliary functions
#--------------------

def create_emission_reduction_dict(path_reduction_txt):
    # read emission reductions per precursor and macro sector
    f = open(path_reduction_txt, 'r')
    emission_reduction_dict = {}
    f.readline()
    while True:
        line = f.readline().rstrip()
        if len(line) == 0:
            break
        value_lst = line.split('\t')
        # sub dictionary per precursor
        precursor = value_lst[0]
        emission_reduction_dict[precursor] = {}
        for snap in range(1, 11):
            emission_reduction_dict[precursor][snap] = float(value_lst[snap]) / 100.0
    f.close()
    
    return emission_reduction_dict

def create_emission_dict(path_emission_cdf, precursor_lst):
    # open the emission netcdf
    rootgrp = Dataset(path_emission_cdf, 'r')
    
    emission_dict = {}
    for precursor in precursor_lst:
        emission_dict[precursor] = rootgrp.variables[precursor][:, :, :]
    
    # get snap, longitude and latitude arrays from emission file
    snap_array = range(1, 11)
    lon_array = rootgrp.variables['longitude'][:]
    lat_array = rootgrp.variables['latitude'][:]
    emission_dict['Nsnaps'] = snap_array
    emission_dict['lon_array'] = lon_array
    emission_dict['lat_array'] = lat_array
    
    # close the emission file
    rootgrp.close()
    
    return emission_dict


# definition of a class managing the low and high resolution emissions aggregation
class Emissions:
    def __init__(self, n_lowres, path_emission_cdf, precursor_lst, path_area_cdf, path_reduction_txt):
        self.path_emission_cdf = path_emission_cdf
        self.path_area_cdf = path_area_cdf
        self.path_reduction_txt = path_reduction_txt
        # create a dictionary with reductions per precursor and macro sector
        emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
        # open the emission netcdf
        emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)
    
        # open the area netcdf
        rootgrp = Dataset(path_area_cdf, 'r')
        reduction_area = rootgrp.variables['AREA'][:] / 100.0
        rootgrp.close()
    
        # calculate a dictionary with the emission reductions per pollutant, macrosector and position
        self.delta_emission_dict = {}
        for precursor in precursor_lst:
            self.delta_emission_dict[precursor] = zeros(emission_dict[precursor].shape)
            # calculate the emission reduction
            # reductions are positive!
            # make the sum over all snap sectors
            for snap in range(1, 11):
                self.delta_emission_dict[precursor][snap - 1, :, :] = emission_dict[precursor][snap - 1, :, :] * reduction_area * emission_reduction_dict[precursor][snap]
        
        # WHEN OUTPUT PER SNAP IS NEEDED, DO NOT SUM YET OVER SNAPS
        
        # sum over all snap sectors
        for precursor in precursor_lst:
            self.delta_emission_dict[precursor] = sum(self.delta_emission_dict[precursor], axis=0)
        
        # self.dem = delta_emission_dict
        self.n_lat = self.delta_emission_dict['PPM'].shape[0]
        self.n_lon = self.delta_emission_dict['PPM'].shape[1]
        self.n_lowres = n_lowres
        self.agg_radius = (n_lowres - 1) / 2
        self.aggregated_emissions = {}  
        self.emissions_hires = {}    
    
    # generates or looks up hi and low resolution emissions for a precursor and coordinate
    def getHiLowEmis(self, precursor, i_lat_target, i_lon_target):
        
        i_lat_mod = i_lat_target % self.n_lowres
        i_lon_mod = i_lon_target % self.n_lowres
        i_tuple_mod = (i_lat_mod, i_lon_mod)        # the modulus of the indexes characterizes the aggregation.
        
        # check if for the modulus of the row-column index tuple some aggregations are done
        if not(i_tuple_mod in self.aggregated_emissions.keys()):
            self.aggregated_emissions[i_tuple_mod] = {}
            self.emissions_hires[i_tuple_mod] = {}
            
            # define the breaks of the aggregation blocks along latitudes and longitudes
            lon_breaks = range((i_lon_target + (self.n_lowres + 1) / 2) % self.n_lowres, self.n_lon, self.n_lowres)
            # add a zero in front and n_lat as last element if necessary
            if lon_breaks[0] > 0:
                lon_breaks = [0] + lon_breaks       # it's a list, not an array!
            if lon_breaks[-1] < self.n_lon:
                lon_breaks = lon_breaks + [self.n_lon]
 
            # define the breaks of the aggregation blocks along latitudes and longitudes
            lat_breaks = range((i_lat_target + (self.n_lowres + 1) / 2) % self.n_lowres, self.n_lat, self.n_lowres)
            # add a zero in front and n_lat as last element if necessary
            if lat_breaks[0] > 0:
                lat_breaks = [0] + lat_breaks       # it's a list, not an array!
            if lat_breaks[-1] < self.n_lat:
                lat_breaks = lat_breaks + [self.n_lat]
            
            self.aggregated_emissions[i_tuple_mod]['n_lat_agg'] = len(lat_breaks) - 1
            self.aggregated_emissions[i_tuple_mod]['n_lon_agg'] = len(lon_breaks) - 1
            self.aggregated_emissions[i_tuple_mod]['lat_breaks'] = lat_breaks
            self.aggregated_emissions[i_tuple_mod]['lon_breaks'] = lon_breaks
            
            # for debugging
#             print(i_tuple_mod)
#             print(lat_breaks)
#             print(lon_breaks)
#             print('\n')
            
        if not(precursor in (self.aggregated_emissions[i_tuple_mod]).keys()):
            n_lat_agg = self.aggregated_emissions[i_tuple_mod]['n_lat_agg']
            n_lon_agg = self.aggregated_emissions[i_tuple_mod]['n_lon_agg']
            lat_breaks = self.aggregated_emissions[i_tuple_mod]['lat_breaks']
            lon_breaks = self.aggregated_emissions[i_tuple_mod]['lon_breaks']
            
            self.aggregated_emissions[i_tuple_mod][precursor] = zeros((n_lat_agg, n_lon_agg))
            self.emissions_hires[i_tuple_mod][precursor] = zeros((self.n_lat, self.n_lon))
            for i_lat_agg in range(n_lat_agg):
                for i_lon_agg in range(n_lon_agg):
                    # select the block to be aggregated from delta_emissions[precursor]
                    emis2aggregate = self.delta_emission_dict[precursor][lat_breaks[i_lat_agg]:lat_breaks[i_lat_agg + 1], lon_breaks[i_lon_agg]:lon_breaks[i_lon_agg + 1]]
                    sum_emissions = emis2aggregate.sum()
                    n_agg = emis2aggregate.size
                    self.aggregated_emissions[i_tuple_mod][precursor][i_lat_agg, i_lon_agg] = sum_emissions
                    # high resolution delta emissions are delta emissions minus average aggregated emissions
                    self.emissions_hires[i_tuple_mod][precursor][lat_breaks[i_lat_agg]:lat_breaks[i_lat_agg + 1], lon_breaks[i_lon_agg]:lon_breaks[i_lon_agg + 1]] = emis2aggregate - sum_emissions / n_agg
#                     print(self.aggregated_emissions[i_tuple_mod][precursor])
#                     print(self.emissions_hires[i_tuple_mod][precursor])
        
        res_dict = {'emis_lowres': self.aggregated_emissions[i_tuple_mod][precursor], 'emis_hires': self.emissions_hires[i_tuple_mod][precursor], \
                    'n_lat_agg': self.aggregated_emissions[i_tuple_mod]['n_lat_agg'], 'n_lon_agg': self.aggregated_emissions[i_tuple_mod]['n_lon_agg']}
        
        return res_dict
                        
# Window class that returns aggregated weighting windows for a given omega
class OmegaWindows:
    def __init__(self, n_lowres, hires_window_size, lowres_window_size):
        self.n_lowres = n_lowres            # number of rows/columns that are aggregated
        self.hires_window_size = hires_window_size
        self.hires_window_radius = (hires_window_size - 1) / 2
        self.hires_inverse_distance = zeros((self.hires_window_size, self.hires_window_size))
        for iw in range(self.hires_window_size):
            for jw in range(self.hires_window_size):
                cell_dist = sqrt((float(iw - self.hires_window_radius)) ** 2 + (float(jw - self.hires_window_radius)) ** 2) 
                self.hires_inverse_distance[iw, jw] = 1 / (1 + cell_dist)  
        
        self.lowres_window_size = lowres_window_size
        self.lowres_window_radius = (lowres_window_size - 1) / 2
        self.lowres_inverse_distance = zeros((self.lowres_window_size, self.lowres_window_size))
        for iw in range(self.lowres_window_size):
            for jw in range(self.lowres_window_size):
                cell_dist = sqrt((float(iw - self.lowres_window_radius)) ** 2 + (float(jw - self.lowres_window_radius)) ** 2) 
                self.lowres_inverse_distance[iw, jw] = 1 / (1 + cell_dist)  
        
        self.aggreg_lowres_window_size = self.lowres_window_size / self.n_lowres
        
        self.hires_omega_windows = {} 
        self.lowres_omega_windows = {}

    def getOmegaWindow(self, omega):
        if omega in self.hires_omega_windows.keys():
            # for this omega the weighting window has been calculated
            pass
        else:
            hires_omega_window = power(self.hires_inverse_distance, omega)
            self.hires_omega_windows[omega] = hires_omega_window
            
            # aggregate the omega window in blocks of n_lowres x n_lowres cells
            lowres_omega_window = zeros((self.aggreg_lowres_window_size , self.aggreg_lowres_window_size))
            for iw in range(self.aggreg_lowres_window_size):
                for jw in range(self.aggreg_lowres_window_size):
                    lowres_omega_window[iw, jw] = (power(self.lowres_inverse_distance[(iw * self.n_lowres):((iw + 1) * self.n_lowres), (jw * self.n_lowres):((jw + 1) * self.n_lowres)], omega)).mean()
            self.lowres_omega_windows[omega] = lowres_omega_window
            
        res_dict = {'omega_window_hires': self.hires_omega_windows[omega], 'omega_window_lowres': self.lowres_omega_windows[omega]} 
        
        return res_dict 
            
# calculate delta concentration for a cell
class srm:
    def __init__(self, path_model_cdf, path_emission_cdf, path_area_cdf, path_reduction_txt):
        # model parameters
        self.n_lowres = 3           # cell size of the side of the square for emission aggregation
        self.hires_window_size = 123  # cell size of the high resolution window (in high resolution cells)
        self.hires_window_radius = (self.hires_window_size - 1) / 2
        self.lowres_window_size = 501 # cell size of the low resolution window (in high resolution cells)
        self.lowres_window_radius = (self.lowres_window_size - 1) / 2
        self.lowres_agg_size = int(self.lowres_window_size / self.n_lowres)
        self.lowres_agg_radius = (self.lowres_agg_size - 1) / 2
        
        # read the model netcdf
        rootgrp = Dataset(path_model_cdf, 'r')
        self.longitude_array = rootgrp.variables['lon'][0, :]
        self.latitude_array = rootgrp.variables['lat'][:, 0]
        self.n_lon = len(self.longitude_array)  # len(rootgrp.dimensions['longitude'])
        self.n_lat = len(self.latitude_array)  # len(rootgrp.dimensions['latitude'])  
        self.inner_radius = int(getattr(rootgrp, 'Radius of influence'))
        self.precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
        alpha_array = rootgrp.variables['alpha'][:, :, :]    
        omega_array = rootgrp.variables['omega'][:, :, :] 
        flatWeight_array = rootgrp.variables['flatWeight'][:, :, :]
        rootgrp.close()     # close model netcdf
    
        # put alpha and omega in a dictionary
        self.alpha_dict = {}
        self.omega_dict = {}
        self.flatWeight_dict = {}
        for i in range(len(self.precursor_lst)):
            self.alpha_dict[self.precursor_lst[i]] = alpha_array[i, :, :]
            self.omega_dict[self.precursor_lst[i]] = omega_array[i, :, :]
            self.flatWeight_dict[self.precursor_lst[i]] = flatWeight_array[i, :, :]
        
        # create a default delta emissions object
        self.delta_emissions = {}
        for precursor in self.precursor_lst:
            self.delta_emissions[precursor] = zeros((self.n_lat, self.n_lon))
        
        # create an emissions class inside the model
        self.emissions = Emissions(self.n_lowres, path_emission_cdf, self.precursor_lst, path_area_cdf, path_reduction_txt)
        
        # create a weigthing window class inside the model
        self.window = OmegaWindows(self.n_lowres, self.hires_window_size, self.lowres_window_size)
        
    # change n_lowres
    def setAggregationParameters(self, n_lowres_new, hires_window_size_new, lowres_window_size_new):
        self.n_lowres = n_lowres_new
        self.hires_window_size = hires_window_size_new 
        self.hires_window_radius = (self.hires_window_size - 1) / 2 
        self.lowres_window_size = lowres_window_size_new 
        self.lowres_window_radius = (self.lowres_window_size - 1) / 2
        
        # update the Emissions class of the model
        self.emissions = Emissions(self.n_lowres, self.emissions.path_emission_cdf, self.precursor_lst, self.emissions.path_area_cdf, self.emissions.path_reduction_txt)
        
        # update the OmegaWindows class
        self.window = OmegaWindows(self.n_lowres, self.hires_window_size, self.lowres_window_size)
                              
    def writeDeltaEmissions(self, path_result_cdf, write_netcdf_output):
        pass
    
    # function that returns delat concentration for a given point    
    def calcDeltaConcCell(self, i_lat, i_lon):
        
        # loop over all precursors
        for precursor in self.precursor_lst:
            # apply averaging window
            alpha_ij = self.alpha_dict[precursor][i_lat, i_lon]
            omega_ij = self.omega_dict[precursor][i_lat, i_lon]
            
            if not(isnan(alpha_ij)):
                # get the weighting windows at high and low resolution
                windows_dict = self.window.getOmegaWindow(omega_ij)

                # get delta emissions at high and low resolution
                emissions_dict = self.emissions.getHiLowEmis(precursor, i_lat, i_lon)
                
                # there are 2 possibilities: padding emissions with zeros or adjusting the selection
                # second option: adjusting selection
                
                # calculate delta_concentration contribution of the high resolution
                # north, south, east and west boundaries of the emissions involved at high resolution.
                north = min(self.n_lat, i_lat + self.hires_window_radius + 1)
                south = max(0, i_lat - self.hires_window_radius)
                west = max(0, i_lon - self.hires_window_radius)
                east = min(self.n_lon, i_lon + self.hires_window_radius + 1)
                
                # north, south, east and west boundaries of the weighting window involved.
                win_east = min(self.hires_window_size, self.hires_window_size - (i_lon + self.hires_window_radius - self.n_lon + 1))
                win_north = min(self.hires_window_size, self.hires_window_size - (i_lat + self.hires_window_radius - self.n_lat + 1))
                win_west = max(0, self.hires_window_radius - i_lon)
                win_south = max(0, self.hires_window_radius - i_lat)
                delta_conc_hires = (emissions_dict['emis_hires'][south:north, west:east] * windows_dict['omega_window_hires'][win_south:win_north, win_west:win_east]).sum() 
                
                
                # calculate delta_concentration contribution of the low resolution
                i_lon_agg = floor((i_lon + 1) / self.n_lowres)    # longitude index of i_lon in low resolution emission field
                i_lat_agg = floor((i_lat + 1) / self.n_lowres)    # latitude index of i_lat in low resolution emission field
                n_lon_agg = emissions_dict['n_lon_agg']
                n_lat_agg = emissions_dict['n_lat_agg']
                
                # north, south, east and west boundaries of the emissions involved at high resolution.
                north = min(n_lat_agg, i_lat_agg + self.lowres_agg_radius + 1)
                south = max(0, i_lat_agg - self.lowres_agg_radius)
                west = max(0, i_lon_agg - self.lowres_agg_radius)
                east = min(n_lon_agg, i_lon_agg + self.lowres_agg_radius + 1)

                # north, south, east and west boundaries of the weighting window involved.
                win_east = min(self.lowres_agg_size, self.lowres_agg_size - (i_lon_agg + self.lowres_agg_radius - n_lon_agg + 1))
                win_north = min(self.lowres_agg_size, self.lowres_agg_size - (i_lat_agg + self.lowres_agg_radius - n_lat_agg + 1))
                win_west = max(0, self.lowres_agg_radius - i_lon_agg)
                win_south = max(0, self.lowres_agg_radius - i_lat_agg)
                delta_conc_lowres = (emissions_dict['emis_lowres'][south:north, west:east] * windows_dict['omega_window_lowres'][win_south:win_north, win_west:win_east]).sum() 
                
                delta_conc = delta_conc_lowres + delta_conc_hires #   + 
                
                ### ADD NO2 MODULE ####
                
                return delta_conc

    # function that returns delat concentration for a given point    
    def calcDeltaConc(self, delta_conc_path):
        delta_conc = zeros((self.n_lat, self.n_lon))
        for i_lat in range(self.n_lat):
            for i_lon in range(self.n_lon):   
                delta_conc[i_lat, i_lon] = self.calcDeltaConcCell(i_lat, i_lon)
        
        # create a netcdf file
        rootgrp = Dataset(delta_conc_path, 'w', format='NETCDF3_CLASSIC')
        
        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', self.n_lat)
        rootgrp.createDimension('longitude', self.n_lon)
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        latitudes[:] = self.latitude_array
        longitudes[:] = self.longitude_array
    
        # create delta concentration data
        delta_conc_pol = rootgrp.createVariable('delta_concentration', 'f4', ('latitude', 'longitude',))
        delta_conc_pol.units = 'ug/m3'
        delta_conc_pol[:] = delta_conc
        
        rootgrp.close()
        
        return delta_conc
        
        


# test function

if __name__ == '__main__':
    
    # input
    input_path = 'D:/workspace/sherpa/trunk/input/20151116_SR_no2_pm10_pm25/'
    path_model_cdf = input_path + 'SR_PM25_Y.nc'
    path_emission_cdf = input_path + 'BC_emi_pm25_Y.nc'
    # path_emission_cdf = 'D:/workspace/sherpa/bug_correction_20161205/one_emission_source.nc'
    path_area_cdf = 'D:/workspace/sherpa/trunk/input/London_region.nc'
    path_reduction_txt = 'D:/workspace/sherpa/trunk/input/user_reduction_snap7.txt'
    
    mySRM = srm(path_model_cdf, path_emission_cdf, path_area_cdf, path_reduction_txt)
    
    start = time()
    delta_conc = mySRM.calcDeltaConc('D:/workspace/sherpa/branches/Accelerator/output/delta_conc_London.nc')
    stop = time()
    print('Calculation time: %f seconds' % (stop - start))
    
    
#     # test window function
#     testwin = OmegaWindows(3, 3, 9)
#     # print(testwin.hires_inverse_distance)
#     savetxt("C:/temp/hires_inverse_distance.csv", testwin.hires_inverse_distance, delimiter="\t")
#     savetxt("C:/temp/lowres_inverse_distance.csv", testwin.lowres_inverse_distance, delimiter="\t")
#     # print(testwin.lowres_inverse_distance)
#     omegawindows = testwin.getOmegaWindow(1)
#     print(omegawindows['hires_omega_window'])
#     savetxt("C:/temp/hires_omega_window.csv", omegawindows['hires_omega_window'], delimiter="\t")
#     savetxt("C:/temp/lowres_omega_window.csv", omegawindows['lowres_omega_window'], delimiter="\t")
#     print(omegawindows['lowres_omega_window'])
#     
#     
#     # test the Emissions class
#     delta_emission_dict = {}
#     delta_emission_dict['PPM'] = ones((380,480))
#     delta_emission_dict['PPM'][5,5] = 10
#     delta_emission_dict['PPM'][0,1] = 10
#     print(delta_emission_dict['PPM'])
#     print(delta_emission_dict['PPM'].sum())
#     
#     test = Emissions(3, delta_emission_dict)
#     start = time()
#     HiLowEmis_dict = test.getHiLowEmis('PPM', 5, 5)
#     HiEmis = HiLowEmis_dict['emis_hires']
#     LowEmis = HiLowEmis_dict['emis_lowres']
#     print(HiEmis.sum())
#     print(LowEmis.sum())
#     stop = time()
#     print('First call: %f' % (stop - start))
#     
#     start = time()
#     emis_low_res = test.getLowResEmis('PPM', 5, 5)
#     stop = time()
#     print('Second call: %f' % (stop - start))
#     
#     print(emis_low_res)
#     emis_low_res = test.getLowResEmis('PPM', 1, 1)
#     print(emis_low_res)