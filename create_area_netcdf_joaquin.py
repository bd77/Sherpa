# =============================================================================
# Created on 25 January 2018
# 
# Create area netcdfs for the source areas used in the JOAQUIN-project.
# A description of the areas can be found in 
# The query 'grid intersect between joaquin and emep or chimere.sql' calculates 
# the fraction of each chimere or emep grid cell in each source area. The result 
# of this query is 'grid intersect between joaquin and emep or chimere.txt'
# 
# =============================================================================


# imports
from netCDF4 import Dataset, stringtochar
from numpy import zeros, ones, array, linspace

# input for the chimere and emep gridding
chimere_dict = {'path_result_cdf': 'D:/SHERPA/JOAQUIN validation/area_netcdfs/chimere/',
                'n_lon': 384 , 'lon_min': -10.4375, 'lon_max': 37.4375,
                'n_lat': 448, 'lat_min': 34.03125, 'lat_max': 61.96875}

emep_dict = {'path_result_cdf': 'D:/SHERPA/JOAQUIN validation/area_netcdfs/emep/',
             'n_lon': 1200 , 'lon_min': -29.95, 'lon_max': 89.95,
             'n_lat': 520, 'lat_min': 30.05, 'lat_max': 81.95}

config_dict = {'chimere': chimere_dict, 'emep': emep_dict}

# path of the grid intersect file
sd = 'D:/SHERPA/JOAQUIN validation/'
grid_intersect_file = sd + 'grid intersect between joaquin and emep or chimere.csv'

# function writing a netcdf with cell percentages (values in [0, 100] belonging to an area
def write_area_netcdf(path_area_cdf, ctm_dict, area_name, area_array):
    # open a new file
    rootgrp = Dataset(path_area_cdf, 'w', format = 'NETCDF3_CLASSIC')
    # create dimensions for latitude, longitude, the cell percentages and the area name
    rootgrp.createDimension('latitude', ctm_dict['n_lat'])
    rootgrp.createDimension('longitude', ctm_dict['n_lon'])
    rootgrp.createDimension('nuts_id', 1)
    rootgrp.createDimension("z", len(area_name))
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    latitude_array = linspace(ctm_dict['lat_min'], ctm_dict['lat_max'], ctm_dict['n_lat'])
    latitudes[:] = latitude_array
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    longitude_array = linspace(ctm_dict['lon_min'], ctm_dict['lon_max'], ctm_dict['n_lon'])
    longitudes[:] = longitude_array
    area = rootgrp.createVariable('AREA', 'f4', ('nuts_id', 'latitude', 'longitude',))
    area[:] = area_array
    NUTS = rootgrp.createVariable('NUTS', 'S1', ('nuts_id', 'z',))
    area_name_out = stringtochar(array([area_name], 'S%d' % (len(area_name))))
    NUTS[0,:] = area_name_out
    rootgrp.close()
    print(path_area_cdf + ' created.')

# read the grid intersect file and put store it in a dictionary
grid_intersect_dict = {}
gif = open(grid_intersect_file, 'r')
header = gif.readline()
while True:
    line = gif.readline().rstrip()
    if len(line) == 0:
        gif.close()
        break
    [ctm, nuts_id, x_cell, y_cell, fraction] = line.split(';')
    # define a cell key lon$lat
    cell_key = x_cell + '$' + y_cell
    if not(ctm in grid_intersect_dict.keys()):
        grid_intersect_dict[ctm] = {}
    if not(nuts_id in grid_intersect_dict[ctm].keys()):
        grid_intersect_dict[ctm][nuts_id] = {}
    grid_intersect_dict[ctm][nuts_id][cell_key] = float(fraction)
            
# print(grid_intersect_dict['chimere']['BE2'])    


for ctm in grid_intersect_dict.keys():
    ctm_dict = config_dict[ctm]
    # output folder for area netcdfs
    path_result_cdf = ctm_dict['path_result_cdf']
    latitude_array = linspace(ctm_dict['lat_min'], ctm_dict['lat_max'], ctm_dict['n_lat'])
    longitude_array = linspace(ctm_dict['lon_min'], ctm_dict['lon_max'], ctm_dict['n_lon'])    
    # rest area
    rest_area_array = ones((1, ctm_dict['n_lat'], ctm_dict['n_lon'])) * 100
    
    # loop over all areas
    for area_name in grid_intersect_dict[ctm].keys():
        # file name for the area netcdf
        path_area_cdf = path_result_cdf + area_name + '.nc'
        # create array with zeros for percentages of cells inside area_name
        area_array = zeros((1, ctm_dict['n_lat'], ctm_dict['n_lon']))
        
        # fill in the fractions
        for cell_key in grid_intersect_dict[ctm][area_name]:
            # get cell latitude and longitude from the key
            [cell_lon, cell_lat] = cell_key.split('$')
            cell_lon = float(cell_lon)
            cell_lat = float(cell_lat)
            pct_in_cell = grid_intersect_dict[ctm][area_name][cell_key] * 100
            
            # get row index of latitude and col index of longitude
            i_lat = 0
            lat_error = float('inf')
            for i in range(len(latitude_array)):
                lat_dist = abs(cell_lat - latitude_array[i])
                if lat_dist < lat_error:
                    lat_error = lat_dist
                    i_lat = i
            # get row index of latitude and col index of longitude
            i_lon = 0
            lon_error = float('inf')
            for i in range(len(longitude_array)):
                lon_dist = abs(cell_lon - longitude_array[i])
                if lon_dist < lon_error:
                    lon_error = lon_dist
                    i_lon = i
            
            area_array[0, i_lat, i_lon] = pct_in_cell 
            # print("lon: %f (%d), lat: %f (%d), pct: %f" % (cell_lon, i_lon, cell_lat, i_lat, pct_in_cell))
            
        # substract area from rest area
        rest_area_array = rest_area_array - area_array
        
        # write the netcdf for the area with 'area_name'
        write_area_netcdf(path_area_cdf, ctm_dict, area_name, area_array)

    # write a netcdf file for the rest area
    path_rest_area_cdf = path_result_cdf + 'joaquin_rest.nc'
    write_area_netcdf(path_rest_area_cdf, ctm_dict, 'joaquin_rest', rest_area_array)