# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:32:03 2018
This module is used to do all the aggregations for the GUI, as a postcompute. 
it produces equivalent files to the old fortran code (but emissions are in Mg)
NB: emissions are a special case. 

STILL MISSING: population weighted: 
@author: peduzem
"""
import pandas as pd
from time import time
from sherpa_auxiliaries import read_nuts_area 
from sherpa_auxiliaries import read_nc

def aggregationbyarea(dct, grd_int_txt, out_path):
    '''Function that aggregates results by area (for both nuts level and fuas)
    
    Inputs: 
        dct: dictionary with the information 
        grd_int_txt: path to the grid intersect file
        out_path: path of the output directory
    Outpus:
        txt files that contain the value
        'NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3' 
        
    '''
    start = time() 
    dct=dct_e
    # read the grid intersect (fua or nuts)    
    area=read_nuts_area(grd_int_txt, calcall=True)
    # levels (they are colled the same in fua and nuts)
    nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']
    
    for nuts_lv in nuts_lvs:
        # prepare pd dataframe for results (index are the names of the 
        # geographical units in each level and keys are the delta and the base
        # case values)
        res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                           columns=dct.keys())
        # for the delta and base case value
        for key in dct.keys():
            # read the corresonding nc
            nc=read_nc(dct[key]['path'])
            # create a pd dataframe which has as columns the delta or base 
            # case value
            nct=pd.DataFrame(columns=[dct[key]['var']])
            # assign to the dataframe the tansposed matrix of
            # the correscponding value 
            nct[dct[key]['var']]=nc.loc[dct[key]['var']].transpose()
            # define the geographical units over which to iterate
            if dct[key]['aggregation'] == 'sume' and nuts_lv == 'NUTS_Lv3':
                # if it is delta emissions, for nuts_lv3 iterate only over
                # the selected areas. 
                arealist=list(pd.read_table('D:/programs/sherpa/app/data/temp/nuts3_selection.txt', 
                                     header=None, sep='\n')[0])
            else: 
                arealist = area[nuts_lv]['area'].index.levels[0]

            for areait in arealist:
            # create a pd dataframe with the areas of each cell 
            # in the geographical unit
                area_a = pd.DataFrame(area[nuts_lv]['area'].loc[areait])
#                area_p = pd.DataFrame(area[nuts_lv]['parea'].loc[areait])
                  
                # if the aggregation mode is average (e.g. concentration levels)
                if dct[key]['aggregation']=='avg':
                    # area averaged values 
                    area_a['mult'] = area_a['area'].multiply(nct.reindex(area_a.index)[dct[key]['var']])
                    areaxvar=area_a[area_a['mult'].notnull()]['mult'].sum()                                      
                    areatot=area_a[area_a['mult'].notnull()]['area'].sum()
                    # if the aggregation mode is sum (e.g. mortality)
                elif dct[key]['aggregation']=='sum':
                    # sum the values in each cell considering the fraction of the 
                    # cell belonging to a geographical entity
                    area_a['mult'] = area_a['area'].multiply(nct.reindex(area_a.index)[dct[key]['var']])
                    areaxvar=area_a[area_a['mult'].notnull()]['mult'].sum()
                    areatot = 1
                elif dct[key]['aggregation']=='sume': 
                    if nuts_lv == 'NUTS_Lv3':
                        area_a = pd.DataFrame(area['ALL_NUTS_Lv0']['area'].loc['ALL_NUTS_Lv0'].loc[area_a.index])
                    area_a['mult'] = area_a['area'].multiply(nct.reindex(area_a.index)[dct[key]['var']])
                    areaxvar=area_a[area_a['mult'].notnull()]['mult'].sum()
                    areatot=1

        
                # this if statement is to take care of areas in the 
                # shape file which are outside the domain (e.g. Reunion)
                if areatot is not 0: 
                    value=areaxvar/areatot
                else:
                    value= float('NaN')
                res[key].loc[areait]=value
        
        
        print('Saving results')
        res['per']=(res['delta'])/res['bc']*100
        res['value']=res['bc']-res['delta']
        res[['value', 'delta', 'per']].to_csv(out_path+nuts_lv, header=True, index=True, sep='\t', na_rep='NaN', mode='w')  
    end = time()
    print('Calculation time  for aggregation', end-start)
    
if __name__ == '__main__': 
    
    out_path='D:/programs/sherpa/app/data/temp/'
    delta_nc='D:/programs/sherpa/app/data/temp/delta_concentration.nc'
    delta_e_nc='D:/programs/sherpa/app/data/temp/delta_emission.nc'
    value_e_nc='D:/programs/sherpa/app/data/temp/delta_emission.nc'
    health_nc='D:/programs/sherpa/app/data/temp/healthimp.nc'
    grd_int_txt='D:/programs/sherpa/app/data/input/models/chimere_7km_fua/selection/grid_intersect'
    invalue_nc='D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_concentrations/BC_conc_PM25_Y.nc'
    invalue_e_nc='D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_emissions/BC_emi_PM25_Y.nc'

    dct_c = {'delta': {'path': delta_nc, 'var':'delta_concentration', 'aggregation':'avg' },
             'bc': {'path': invalue_nc, 'var':'conc', 'aggregation':'avg' }}

    dct_e = {'delta': {'path': delta_e_nc, 'var':('PPM','Nsnaps07'), 'aggregation':'sume' },
             'bc': {'path': invalue_e_nc, 'var':('PPM','Nsnaps07'), 'aggregation':'sum' }}
    dct_hi = {'delta': {'path': health_nc, 'var':'d_mort', 'aggregation':'sum'},
              'bc': {'path': health_nc, 'var':'v_mort', 'aggregation':'sum'}}
#    
    aggregationbyarea(dct_e, grd_int_txt, out_path)