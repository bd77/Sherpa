'''
Created on 16 Jun 2017

create area netcdfs for the 112 fuas and the Chimere grid

@author: degraba
'''

import psycopg2
from netCDF4 import Dataset, stringtochar
from numpy import zeros, array
from os import path

# imput for the chimere gridding
chimere_dict = {'path_result_cdf': 'D:/SHERPA/FUA112/fua_area_cdfs/allAreas_chimere/',
                'example_area_netcdf': 'D:/SHERPA/FUA112/fua_area_cdfs/allAreas_chimere/Amsterdam_City.nc'}

emep_dict = {'path_result_cdf': 'D:/SHERPA/FUA112/fua_area_cdfs_emep/',
             'example_area_netcdf': 'O:/Integrated_assessment/SHERPA/20170622_emep_first_results/1_base_emissions/BC_emi_PM25_Y.nc'}

config_dict = {'chimere': chimere_dict, 'emep': emep_dict}

for grid_name in ['chimere']: # config_dict.keys():
    
    # output folder for area netcdfs
    path_result_cdf = config_dict[grid_name]['path_result_cdf']
    
    # database connection
    host_name = "localhost"
    database_name = "sherpa"
    user_name = "postgres"
    pw = "root"
    schema = 'sherpa'
    tbl_fua = 'fua112'
    tbl_grid = 'grids'
    
    try:
        conn = psycopg2.connect(host=host_name, database=database_name, user=user_name, password=pw);
    except:
        print("I am unable to connect to the database")
    cur = conn.cursor()
    
    # read an example netcdf
    example_area_netcdf = config_dict[grid_name]['example_area_netcdf']
    rootgrp_example = Dataset(example_area_netcdf, 'r')
    longitude_array = rootgrp_example.variables['longitude'][:]
    latitude_array = rootgrp_example.variables['latitude'][:]
    # n_nuts = len(rootgrp_example.dimensions['nuts_id'])
    n_lon = len(longitude_array)  
    n_lat = len(latitude_array)  
    # nuts_codes_raw = rootgrp_example.variables['NUTS'][:]
    # nuts_codes = []
    rootgrp_example.close()
    
    # get a list of all areas: 112 cities, their commuting zones and background per country
    query = "SELECT id, codeid FROM sherpa.fua112"
    cur.execute(query)
    rows = cur.fetchall()
    area_dict = {}
    for row in rows:
        area_id = row[0]
        area_name = row[1]
        area_dict[area_id] = area_name
    
    
    query = """SELECT a.name_asci, st_x(st_centroid(g.the_cell)) AS x, st_y(st_centroid(g.the_cell)) AS y, 
                st_area(st_intersection(a.geom, g.the_cell))/st_area(g.the_cell) AS pct_in_cell
                FROM (SELECT * FROM sherpa.fua112 WHERE id = %d) a, 
                     (SELECT g.* 
                      FROM sherpa.grids g, 
                           (SELECT ST_SetSRID(st_extent(a.geom), 4326) AS box FROM sherpa.fua112 a WHERE id = %d) b
                           WHERE st_intersects(g.the_cell, b.box) = true AND g.grid_name = '%s') g;
            """
    
    for area_id in area_dict.keys():
        # name of the netcdf file
        path_area_cdf = path_result_cdf + area_dict[area_id] + '.nc'
        path_area_cdf = path_area_cdf.replace('\r', '') # only for Atina_city, it contains a carriage return
        area_name = area_dict[area_id]
        
        # if it already exists skip and go to next
        if path.exists(path_area_cdf):
            print(path_area_cdf + ' exists')
        else:
            print('Creating ' + path_area_cdf)
            # substitute the fua id in the query
            fua_query = query % (area_id, area_id, grid_name)
            print(fua_query)
            # execute the query: get the share of the area of each cell inside the source area (for those cells inside the boundary box of the source area
            cur.execute(fua_query)
            rows = cur.fetchall()
            
            # create an area netcdf for the source area
            rootgrp = Dataset(path_area_cdf, 'w', format = 'NETCDF3_CLASSIC')
            rootgrp.createDimension('latitude', n_lat)
            rootgrp.createDimension('longitude', n_lon)
            rootgrp.createDimension('nuts_id', 1)
            rootgrp.createDimension("z", len(area_name))
            latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
            latitudes.units = "degrees_north"
            latitudes[:] = latitude_array
            longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
            longitudes.units = "degrees_east"
            longitudes[:] = longitude_array
            area = rootgrp.createVariable('AREA', 'f4', ('nuts_id', 'latitude', 'longitude',))
            NUTS = rootgrp.createVariable('NUTS', 'S1', ('nuts_id', 'z',))
            area_name_out = stringtochar(array([area_name], 'S%d' % (len(area_name))))
            NUTS[0,:] = area_name_out
            array_pct_in_area = zeros((1, n_lat, n_lon))
            
            
            for row in rows:
                fua_name = row[0]
                target_cell_lon = row[1]
                target_cell_lat = row[2]
                pct_in_cell = row[3]
                
                if pct_in_cell > 0:
                    # get row index of latitude and col index of longitude
                    i_lat_target = 0
                    lat_error = float('inf')
                    for i in range(len(latitude_array)):
                        lat_dist = abs(target_cell_lat - latitude_array[i])
                        if lat_dist < lat_error:
                            lat_error = lat_dist
                            i_lat_target = i
                    
                    i_lon_target = 0
                    lon_error = float('inf')
                    for i in range(len(longitude_array)):
                        lon_dist = abs(target_cell_lon - longitude_array[i])
                        if lon_dist < lon_error:
                            lon_error = lon_dist
                            i_lon_target = i
            
                    array_pct_in_area[0, i_lat_target, i_lon_target] = pct_in_cell * 100
                    
                    # print("lon: %f, lat: %f, pct: %f" % (target_cell_lon, target_cell_lat, pct_in_cell))
            
            area[:] = array_pct_in_area  
            NUTS = fua_name  
            rootgrp.close()
    


if __name__ == '__main__':
    
    pass