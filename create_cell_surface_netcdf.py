'''
Created on 12 Jul 2017

Create a netcdf with the cell surface for the emep and chimere grid.
These surfaces are needed to calculate absolute emissions.
WARNING; emep and chimere use different units!!!

@author: degraba
'''

from netCDF4 import Dataset

input_dict = {'emep': {'surface_input_file': 'O:/Integrated_assessment/EMEP-10km/input/gridarea/GridArea_m2.nc',
                       'surface_varname': 'GridArea_m2'},
              'chimere': {'surface_input_file': 'O:/Integrated_assessment/CHIMERE-7km/areas and emissions/JRC01.nc',
                          'surface_varname': 'surface'}}
output_path = 'D:/SHERPA/FUA112/cell_surface_cdfs/'

for model in input_dict.keys():

    # open netcdf with cell surfaces, read the surface array
    surface_cdf = Dataset(input_dict[model]['surface_input_file'], 'r')
    longitude_array = surface_cdf.variables['lon'][:]
    latitude_array = surface_cdf.variables['lat'][:]
    n_lon = len(longitude_array)  
    n_lat = len(latitude_array)  
    surface_array = surface_cdf.variables[input_dict[model]['surface_varname']][:]
    surface_units = surface_cdf.variables[input_dict[model]['surface_varname']].units
    surface_cdf.close()
    
    # write a netcdf with the cell surfaces, a lat and lon array and units
    surface_cdf_file = output_path + model + '_cell_surface.nc'
    rootgrp = Dataset(surface_cdf_file, 'w', format = 'NETCDF3_CLASSIC')
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    latitudes[:] = latitude_array
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    longitudes[:] = longitude_array
    surface_var_new = rootgrp.createVariable('surface', 'f4', ('latitude', 'longitude',))
    if model == 'chimere':
        surface_var_new[:] = surface_array
        surface_var_new.units = surface_units
    if model == 'emep':
        surface_var_new[:] = surface_array / 1e6    # convert m2 into km2
        surface_var_new.units = 'km2'
            
    rootgrp.close()



if __name__ == '__main__':
    pass