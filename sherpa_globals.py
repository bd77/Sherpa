'''
Created on Jul 14, 2015

define some global variables

@author: degraba
'''

# variabels for testing

# netcdf with cells where reductions have to be applied (value between 0 and 1)
# path_area_cdf_test = 'input/EMI_RED_NUTS3_ITALY.nc'
path_area_cdf_test = 'input/London_region.nc'
# reductions per precursor and macro sector
path_reduction_txt_test = 'input/user_reduction_snap7.txt'
path_reduction50all_txt_test = 'input/user_reduction_all50.txt'
# reductions per precursor and macro sector for module 3a and 3b
path_reduction_mod3a1P_txt_test = 'input/potency_reduction_module3a1P.txt'
path_reduction_mod3a2P_txt_test = 'input/potency_reduction_module3a2P.txt'
path_reduction_mod3b_txt_test = 'input/potency_reduction_module3b.txt'
# netcdf with model parameters per cell
#path_model_cdf_test = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
path_model_cdf_test = 'input\\20170322_v18_SrrResults_PotencyBased\\3_source_receptors\\SR_PM25_Y_20170322_potencyBased.nc'
# folder where output will be put
path_result_cdf_test = 'output/'
# progress log is used when module 1 is called by another module
path_nuts0_cdf_test = 'input/EMI_RED_ATLAS_NUTS0.nc'
path_nuts1_cdf_test = 'input/EMI_RED_ATLAS_NUTS1.nc'
path_nuts2_cdf_test = 'input/EMI_RED_ATLAS_NUTS2.nc'
path_nuts3_cdf_test = 'input/EMI_RED_ATLAS_NUTS3.nc'

# absolute emission per cell and macro sector
#path_emission_cdf_test = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
#path_base_conc_cdf_test = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
#path_dust_conc_cdf_test = 'input\\pDUST-pSALT\\pDUST-25-basecase.nc'
#path_salt_conc_cdf_test = 'input\\pDUST-pSALT\\pSALT-25-basecase.nc'
#
path_emission_cdf_test = 'input\\20170322_v18_SrrResults_PotencyBased\\1_base_emissions\\BC_emi_PM25_Y.nc'
path_base_conc_cdf_test = 'input\\20170322_v18_SrrResults_PotencyBased\\2_base_concentrations\\BC_conc_PM25_Y.nc'
path_base_NO2conc_cdf_test = 'input\\20170322_v18_SrrResults_PotencyBased\\2_base_concentrations\\BC_conc_NO2_NO2eq_Y_mgm3.nc'

path_dust_conc_cdf_test = 'input\\pDUST-pSALT\\pDUST-25-basecase.nc'
path_salt_conc_cdf_test = 'input\\pDUST-pSALT\\pSALT-25-basecase.nc'

path_healthbl_test = 'input/impacts/healthbl_nc.nc'
path_config_json_test = 'config/sharedvariables.json'
fua_intersect_dir = 'input/selection/gridnew/fua/'
nuts_intersect_dir = 'input/selection/gridnew/nuts/'
dbf_dir = 'input/selection/gridnew/'
target_list = 'input/AM_targets.txt'
path_natural_dir_test = 'input/pDUST-pSALT/'
aggr_zones='nuts'

path_pop_mat_test = 'input\\population.mat'
path_grid_txt = 'input\\selection\\grid_intersect'
gcities_txt = 'input\\selection\\grid_int_gcities'
fua_txt = 'input\\selection\\grid_int_fua'
path_mortbaseline = 'input\\impacts\\mortalitybaseline.xlsx'
path_tiff0 = 'input\\pop\\7km_Qi_2010.tif' # population distribution according to LUISA
path_tiff = 'input\\pop\\Qi_2015_mosaic.tif'
path_healthbl = 'input\\impacts\\healthbl_nc.nc'
path_healthbl0 = 'input\\impacts\\healthbl0_nc.nc'
json_path = 'config\\sharedvariables.json'
# list of precursors
# order important, it's the order in the alpha and omega arrays
# precursor_lst = ['NOx', 'NMVOC', 'NH3', 'PM25', 'SOx']

# order important, it's the order in the alpha and omega arrays
sector_lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]   #

# fixed reduction percentage for potency calculation
alpha_potency = float(50)
