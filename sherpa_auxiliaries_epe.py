'''
Created on April 10, 2017

auxiliary functions for SHERPA, to be added to the sherpa_auxiliaries file

@author: peduzem
'''

from netCDF4 import Dataset
import numpy as np  # for scientific operators
from osgeo import gdal, ogr, osr  #conda install -c conda-forge gdal
import pickle  # to save results direclty as python objects
import pandas as pd
from module7_custom import read_nuts_area
from sherpa_globals import path_model_cdf_test, grid_txt, gcities_txt, fua_txt, path_emission_cdf_test

from sherpa_auxiliaries import (create_emission_dict)

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

# -----------------------------------------------------------------------------
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
# -----------------------------------------------------------------------------
def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
# -----------------------------------------------------------------------------

def gridint_toarray(level, parea, code):
    """
    Reads the grid intersect txt files and creates an array with the specified
    dimensions with the fraction of each cell beleonging to the specified area
  . Needs to import the function from DENISEP
    from module7_custom import read_nuts_area

    INPUT:
        level = 'NUTS_Lv0' or GCITY_CODE or FUA_CODE
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

    nuts_info = read_nuts_area(grid_txt, calcall=True)
    nuts_info.update(read_nuts_area(gcities_txt, nullnut='LAND000'))
    nuts_info.update(read_nuts_area(fua_txt, nullnut='LAND000'))

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


def write_nc(array, path_nc, name_var, unit_var, addnutsid=False):
    ''' Function to write an array in a netcdf file,
        input:
            - array: data to write
            - path_nc: path of netcdf file
            - name_var: name for data in array
            - unit_var: units for data in array
            - addnutsid: if True the layer nuts_id is added so that the
                nectcdf file is consistent with the ones provided
                by terraria
    @author: peduzem
    '''
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()
    fh = Dataset(path_nc, mode='w', format='NETCDF3_CLASSIC')
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
        fh.variables[name_var].units = unit_var
    fh.close()


def bm_precompute(perreduction, sources):
    """
    Precompute function for the blame matrix.
    Input:
        - perreduction: percentage reduction of precursor emission
        - sources: list of countries
    Output:
        -dfemidelta: dataframe with the total emission reduction [Mg]
        It also prdouces:
        -reduction_txt (for each precursor)
        -area_nc: for each source/target country

    Created on Thu Apr 20 16:26:42 2017
    @author: peduzem
    """
    # -------------------------------------------------------------------------
    # Write reduction text
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()

    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
    for precursor in precursor_lst:
        df_red = pd.DataFrame(np.zeros(
                              shape=(len(precursor_lst), len(ms_list))),
                              index=precursor_lst, columns=ms_list)
        path_reduction_txt = 'workdir\\redbm_{}.txt'.format(precursor)
        df_red.ix[precursor] = perreduction
        df_red.to_csv(path_reduction_txt, sep='\t', index_label='POLL')
    # -------------------------------------------------------------------------
    level = 'NUTS_Lv0'
    parea = 'parea'
    path_surf_nc = 'input/JRC01.nc'
    # Total emissions by country and precursor
    # -------------------------------------------------------------------------
    # create a dictionary with emissions per precursor, macrosector and postion
    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)
    # Area of each cell in the domain
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    # area_units = rootgrp.variables['surface'].units
    rootgrp.close()
    # dataframe with emissions and correponding emission reduction
    # by country and pollutant
    dfemi = pd.DataFrame(index=sources, columns=precursor_lst)
    for source in sources:
        path_area_nc = 'workdir/area_{}.nc'.format(source)
        area = gridint_toarray(level, parea, source)
        write_nc(area, path_area_nc, 'AREA', '%')
        emi = {precursor: np.sum(
               np.sum((emission_dict[precursor]*surf), axis=0) * area/100)
               for precursor in precursor_lst}
        dfemi.ix[source] = emi
        dfemidelta = dfemi * perreduction / 100
    dfemi.to_csv('workdir\\emi.csv')
    dfemidelta.to_csv('workdir\\emidelta.csv')
    save_obj(dfemidelta, 'dfemidelta')
    return dfemidelta
#if __name__ == '__main__':
#    path_tiff = 'input/pop/7km_Qi_2010.tif'
#    popall = tiftogridgeneral(path_tiff)