from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/  
import numpy as np
import pandas as pd  
'''
NAME
    Reads netcdf files and put results in a dataframe with columns for gridpoints and (evenutally) ultindexed rows
PURPOSE
    Reads netcdf files and put results in a dataframe with columns for gridpoints and (evenutally) ultindexed rows
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
def read_nc_df(nc_file):
    #nc_file='input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    #nc_file='input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    nc_data = Dataset(nc_file, 'r') 
    nc_dims = [dim for dim in nc_data.dimensions]
    nc_vars = [var for var in nc_data.variables]
    #sometimes the latitude is written just with lat as in model data
    latname=list(filter(lambda x: x in nc_vars, ['Lat','lat','latitude']))[0]
    #latname=list(filter (lambda x: 'lat' in x, nc_vars))[0]
    lats = nc_data.variables[latname][:] 
    lonname=list(filter(lambda x: x in nc_vars, ['Lon','lon','longitude']))[0]
    #lonname=list(filter (lambda x: 'lon' in x, nc_vars))[0]
    lons = nc_data.variables[lonname][:] 
    #if there are three dimensional arrays
    if len(nc_dims)==3:
        ncz=str(list(set(nc_dims)-set(['latitude','longitude']))[0])
        nz=range(len(nc_data.dimensions[ncz]))
        if ncz=='pollutant':
            strpoll=nc_data.Order_Pollutant
            nznames=strpoll.split(', ')
        else:
            nznames=[ncz + s for s in map(str,range(1,len(nc_data.dimensions[ncz])+1))]         
    #create an index with lat and lon
    #latrep=map(str, np.repeat(lats,len(lons)))
    #lonrep=map(str, np.tile(lons,len(lats)))
    #trasform variables arrays in vectors
    #allvar={'lat_lon':map(lambda (x,y): x+'_'+y, zip(latrep, lonrep))}
    #create lat and lon info
    if len(lats.shape)==2 and len(lons.shape)==2:
        nrow=lats.shape[0] 
        ncol=lats.shape[1]
        lon=lons.ravel()
        lat=lats.ravel()
    else:
        nrow=len(lats)
        ncol=len(lons)
        lon=np.tile(lons,nrow)
        lat=np.repeat(lats,ncol)
         
    y=np.repeat(range(1, nrow+1),ncol)
    x=np.tile(range(1, ncol+1),nrow) 
    row=list(map(str,y))
    col=list(map(str,x))
    index_grid=list(map(lambda x: '_'.join(x),list(zip(col,row))))

    allvar={}
    allvar['coord']=pd.DataFrame(lon,columns=['lon'])
    allvar['coord']['lat']=lat
    allvar['coord']['x']=x
    allvar['coord']['y']=y
    allvar['coord'].index=index_grid
    nc_vars.remove(latname)
    nc_vars.remove(lonname)
    for var in nc_vars:
        varnc=nc_data.variables[var][:]
        if len(nc_dims)==3:
            allvar[var]=pd.concat(map(lambda sn : pd.Series(varnc[sn].ravel()),nz),axis=1)
            allvar[var].columns=nznames
        else:
             allvar[var]=pd.DataFrame(varnc.ravel())
             allvar[var].columns=[var]
        allvar[var].index=index_grid
        #allvarnc[var]=allvarnc[var].transpose()
        #index_var = pd.MultiIndex.from_tuples(zip(np.repeat(var,len(nz)),nznames), names=['vars', ncz])
        #allvar[var].columns=index_var
    reform = {(outerKey, innerKey): values for outerKey, innerDict in allvar.items() for innerKey, values in innerDict.items()}
    df=pd.DataFrame(reform)  
    return df.transpose()
'''
NAME
    Write a 2D matrix or a dictionary of 2D matrices on a netcdf file
PURPOSE
    CWrite a 2D matrix or a dictionary of 2D matrices on a netcdf file
PROGRAMMER(S)
    Denise Pernigotti starting from Bart routine for writing delta conc netcdf in module 1
REVISION HISTORY
    
REFERENCES
    
''' 
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/     
def write_dict_nc(dc_dic,lats,lons,unit,filenc):
    # create a result netcdf
    varnames=list(dc_dic.columns)
    dfmatnc=dict(zip(varnames,list(map(lambda x: df2mat(dc_dic[x]),varnames))))
    rootgrp = Dataset(filenc+'.nc', 'w', format='NETCDF3_CLASSIC')     
    # create dimensions in the netcdf file
    rootgrp.createDimension('latitude', len(lats))
    rootgrp.createDimension('longitude', len(lons))
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    latitudes[:] = lats
    longitudes[:] = lons
    for var in  varnames: 
        # create delta concentration data
        varnc = rootgrp.createVariable(var, 'f4', ('latitude', 'longitude',))
        varnc.units = unit
        varnc[:] =dfmatnc[var]
        
    rootgrp.close()
    
'''
NAME
    Reads SHERPA inputs and manage them in order to minimize the dataframes/arrays dimensions
PURPOSE
    Input the path to netcdf files, the modelled pollutant name (PM25,PM10,...), the emission_ponts selected in order to exclude points that have not 
    been cosidered in the optimization (e.g. sea ponts), and an option to eventually include natural emissions for PM. The output also inclide a dictionary 
    with 3D array in 'data', and dimensions labels in 'dimnames'
PROGRAMMER(S)
    Denise Pernigotti, March 2017
REVISION HISTORY
    
REFERENCES
    
'''    
def clean_input(nc_path,pollname,emitting_points=[],include_natural=False):
       #file to be imported for SR model
    emission_nc = nc_path + 'base_emissions\\BC_emi_'+pollname+'_Y.nc'
    concentration_nc = nc_path + 'base_concentrations\\BC_conc_'+pollname+'_Y.nc'
    model_nc = nc_path + 'source_receptors\\SR_'+pollname+'_Y.nc'
 
    #read netcdf files, put the data in multindex dataframes and check consistency 
    emissions = read_nc_df(emission_nc) 
    concentration = read_nc_df(concentration_nc) 
    model= read_nc_df(model_nc) 
    if include_natural:
        if pollname[0:2]=='PM':pmsize=pollname[2:4]
        salt_nc = 'input\\pDUST-pSALT\\pSALT-'+pmsize+'-basecase.nc'
        dust_nc = 'input\\pDUST-pSALT\\pDUST-'+pmsize+'-basecase.nc'
        salt= read_nc_df(salt_nc)
        dust= read_nc_df(dust_nc)
        
    #check consistency of model and emissions
    if model.loc['coord'].equals(emissions.loc['coord']): 
        print ('OK latitude and longitude in matrices emissions and model are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")    


    #check consistency of model and concentrations
    if model.loc['coord'].equals(concentration.loc['coord']): 
        print ('OK latitude and longitude in matrices model and conc are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")         
    
    #save coord info
    coordinates=emissions.loc['coord',].transpose()
    lat_vec=np.unique(coordinates['lat'])
    lon_vec=np.unique(coordinates['lon'])
    #find grid resolution in degrees 
    dlon_res=min(lon_vec[1:]-lon_vec[:-1])
    dlat_res=min(lat_vec[1:]-lat_vec[:-1])
    
    #remove coord columns
    model=model.drop('coord',level=0)
    emissions=emissions.drop('coord',level=0)
    concentration=concentration.drop('coord',level=0)
    #fake sum to remove one level
    concentration=concentration.groupby(level=[0]).sum()  
    
    inner_radius=False
    if 'flatWeight' in model.index.get_level_values(0):
        inner_radius = int(getattr(Dataset(model_nc,'r'), 'Radius of influence'))
        
    #remove points with zero emissions
    nonzero_em_idx=emissions.columns[emissions.groupby(level=0).sum().sum()>0]
    emissions=emissions[nonzero_em_idx]
    print('keeping ',len(nonzero_em_idx), ' grid points with non zero emissions')

    #if sea and other areas are not modelled land (both as source and receptor) exclude them
    if len(emitting_points)>0: 
        emissions_idx=list(set(emissions.columns).intersection(emitting_points))
        emissions=emissions[emissions_idx]
        print('keeping ',len(emissions_idx), ' only modelled grid points for emissions')

    #remove not modelled grid points
    model_sum=model.loc['alpha'].sum()
    model_points=model_sum[model_sum>0].index
    print ("in SR data keep only ",len(model_points)," where alpha is non zero")
    model=model[model_points]

    #remove not modelled precursors
    model_sumgrid=model.loc['alpha'].sum(axis=1)
    precursors=list(model_sumgrid[model_sumgrid>0].index)
    print ("in model data keep only ",list(precursors)," modelled precursors")
    model=model.loc[(slice(None),precursors),:]  
    emissions=emissions.loc[(precursors,slice(None)),:]
     
    #snaps names and precursors names from emissions
    snaps=emissions.index.levels[1]
    snaps=list(snaps[~snaps.isin(['lat','lon','x','y'])])

    #Convert the emissions dataframe to a 3D array for calculation speed
    emissions_array={}
    emissions_idx=emissions.columns
    emissions_array['dimnames']={}
    emissions_array['dimnames']['snaps']=snaps
    emissions_array['dimnames']['emissions_idx']=emissions_idx
    emissions_array['dimnames']['precursors']=precursors
    emissions_array['data']=np.zeros(shape=(len(snaps),len(emissions_idx),len(precursors)))
    for isp,prec in enumerate(precursors):
        emissions_array['data'][:,:,isp]=np.array(emissions.loc[(prec,slice(None)),:])  
        
    outdic={'coordinates':coordinates,'dlat_res':dlat_res,'dlon_res':dlon_res,'model':model,'concentration':concentration,
    'emissions':emissions,'emissions_array':emissions_array,'inner_radius':inner_radius}
    return outdic

