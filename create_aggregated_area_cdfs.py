'''
Created on 30 Jun 2017

Create 4 aggregated area netcdfs per city:
- <city>_City (not an aggregation, exist already so not calculated or copied)
- <city>_Comm (not an aggregation, exist already so not calculated or copied)
- <city>_National (Other national emissions except City and Commuting zone, calculated for each city)
- <city>_International named <Country>_International (All international contributions, same for all cities of the same country)

@author: degraba
'''

import os
from netCDF4 import Dataset, stringtochar
from numpy import array

root_folder = 'D:/SHERPA/FUA112/'
inventory_list = ['emep', 'chimere']
area_cdf_folder_format = root_folder + 'fua_area_cdfs_%s/' 
agg_area_cdf_folder_format = root_folder + 'agg_fua_area_cdfs_%s/' 


# read city list
city_list_file = root_folder + 'city_list_fua112.txt'
city_dict = {}
fcity = open(city_list_file, 'r')
fcity.readline()     # read header
city_dict = {}
while True:
    line = fcity.readline().rstrip()
    if len(line) == 0:
        break
    [cityname, lat, lon, codeid, country_code] = line.split(';')
    city_dict[cityname] = {'codeid': codeid, 'lat': float(lat), 'lon': float(lon), 'country_code': country_code}

for inventory in inventory_list:
    print('+++++++++%s+++++++++++++' % (inventory))
    # folder with the netcdfs to be aggregated
    area_cdf_folder = area_cdf_folder_format % (inventory)
    source_area_list = os.listdir(area_cdf_folder)
    source_area_dict = {}
    # remove extension '.nc' from area name
    for i in range(len(source_area_list)):
        source_area_name = source_area_list[i].replace('.nc', '')
        if source_area_name == 'Background_BE':
            print(source_area_name)
        
        source_area_list[i] = source_area_name
        [p1, p2] = source_area_name.split('_')
        if p1 == 'Background':
            source_area_country = p2
        else:
            source_area_city = p1
            source_area_country = city_dict[source_area_city]['country_code']
        source_area_dict[source_area_name] = source_area_country 
    # print(source_area_dict)
    # output folder for aggregated netcdfs
    agg_area_cdf_folder = agg_area_cdf_folder_format % (inventory)
    
    # loop over all cities
    for target_city in city_dict.keys():
        print('target_city: %s' % (target_city))
        target_country = city_dict[target_city]['country_code']
        
        # look for national areas excluding the City and Commuting area
        area_national = 0
        area_international = 0
        for source_area in source_area_dict.keys():
            
            if source_area == target_city + '_City' or source_area == target_city + '_Comm':
                pass # print('%s not included' % (source_area))
            elif source_area_dict[source_area] == target_country:
                # print('%s included in national' % (source_area))
                source_area_netcdf = area_cdf_folder + source_area + '.nc'
                rootgrp_source = Dataset(source_area_netcdf, 'r')
                longitude_array = rootgrp_source.variables['longitude'][:]
                latitude_array = rootgrp_source.variables['latitude'][:]
                n_lon = len(longitude_array)  
                n_lat = len(latitude_array)
                area = rootgrp_source.variables['AREA'][:] 
                rootgrp_source.close()
                area_national = area_national + area
            else: 
                # print('%s included in international' % (source_area))
                source_area_netcdf = area_cdf_folder + source_area + '.nc'
                rootgrp_source = Dataset(source_area_netcdf, 'r')
                longitude_array = rootgrp_source.variables['longitude'][:]
                latitude_array = rootgrp_source.variables['latitude'][:]
                n_lon = len(longitude_array)  
                n_lat = len(latitude_array)
                area = rootgrp_source.variables['AREA'][:] 
                rootgrp_source.close()
                area_international = area_international + area
                
        # write an area netcdf for the national and international areas
        national_area_name = target_city + '_National'
        national_area_cdf = agg_area_cdf_folder + national_area_name + '.nc'
        rootgrp = Dataset(national_area_cdf, 'w', format = 'NETCDF3_CLASSIC')
        rootgrp.createDimension('latitude', n_lat)
        rootgrp.createDimension('longitude', n_lon)
        rootgrp.createDimension('nuts_id', 1)
        rootgrp.createDimension("z", len(national_area_name))
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        latitudes[:] = latitude_array
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        longitudes[:] = longitude_array
        area = rootgrp.createVariable('AREA', 'f4', ('nuts_id', 'latitude', 'longitude',))
        NUTS = rootgrp.createVariable('NUTS', 'S1', ('nuts_id', 'z',))
        area_name_out = stringtochar(array([national_area_name], 'S%d' % (len(national_area_name))))
        NUTS[0,:] = area_name_out
        area[:] = area_national  
        NUTS = national_area_name  
        rootgrp.close()
        
        international_area_name = city_dict[target_city]['country_code'] + '_International'
        if not(os.path.exists(international_area_name)):
            international_area_cdf = agg_area_cdf_folder + international_area_name + '.nc'
            rootgrp = Dataset(international_area_cdf, 'w', format = 'NETCDF3_CLASSIC')
            rootgrp.createDimension('latitude', n_lat)
            rootgrp.createDimension('longitude', n_lon)
            rootgrp.createDimension('nuts_id', 1)
            rootgrp.createDimension("z", len(international_area_name))
            latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
            latitudes.units = "degrees_north"
            latitudes[:] = latitude_array
            longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
            longitudes.units = "degrees_east"
            longitudes[:] = longitude_array
            area = rootgrp.createVariable('AREA', 'f4', ('nuts_id', 'latitude', 'longitude',))
            NUTS = rootgrp.createVariable('NUTS', 'S1', ('nuts_id', 'z',))
            area_name_out = stringtochar(array([international_area_name], 'S%d' % (len(international_area_name))))
            NUTS[0,:] = area_name_out
            area[:] = area_international  
            NUTS = international_area_name  
            rootgrp.close()

            
            

if __name__ == '__main__':
    pass