# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:32:03 2018
This module is used to do all the aggregations for the GUI, as a postcompute. 
it produces equivalent files to the old fortran code (but emissions are in Mg)

NB: emissions are a special case. 

STILL MISSING: average over threshold 
@author: peduzem
"""
import pandas as pd
from time import time
import json
import ast 


#import os as os 

from sherpa_auxiliaries import read_nuts_area 
from sherpa_auxiliaries import read_nc

def aggregationbyarea(aggrinp_txt):
    '''Function that aggregates results by area (for both nuts level and fuas)
    
    Inputs: 
        -aggrinp_txt: txt file with the information:
        {
        "delta": {
                "path": "D:/programs/sherpa/app/data/temp/delta_concentration.nc",
                "var": "delta_concentration",
                "aggregation": "('avg','pop')"},
        "bc": {
             "path": "D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_concentrations/BC_conc_PM25_Y.nc",
             "var": "conc",
             "aggregation": "('avg','pop')"},
     	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
    	"output-dir":"D:/programs/sherpa/app/data/temp/"
        }
    

    Outpus:
        -txt files that contain the value
        'NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3' 
        
    '''
#    aggrinp_txt='D:/programs/sherpa/app/data/temp/aggrinp.txt'
    json_file = open(aggrinp_txt)
    json_str = json_file.read()
    dct = json.loads(json_str)
    grd_int_txt = dct['grid-intersect']
    out_path =  dct['output-dir']
    start = time()    # read the grid intersect (fua or nuts)    
    area=read_nuts_area(grd_int_txt, calcall=True)
    # levels (they are colled the same in fua and nuts)
    nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']
    dct_ms={'Nsnaps01':'MS1','Nsnaps02':'MS2','Nsnaps03':'MS3',
            'Nsnaps04':'MS4','Nsnaps05':'MS5','Nsnaps06':'MS6',
            'Nsnaps07':'MS7','Nsnaps08':'MS8', 'Nsnaps09':'MS9', 
            'Nsnaps10':'MS10'}
    grd_int=pd.DataFrame.from_csv(grd_int_txt+'.txt', sep='\t')
    # if we are aggregating emission reductions we need to consider the 
    # exact area defined by the user or reductions are smeared out on the grid. 
    if ast.literal_eval(dct['bc']['aggregation'])=='sume':
        nuts_lv = 'NUTS_Lv3'
        # get tuple defining precursor and sector
        tpl = ast.literal_eval(dct['bc']['var'])
        # read reduction file 
        red = pd.read_table(out_path+'user_reduction.txt', index_col=['POLL'])
        # define df for results
        res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                           columns=['bc','delta'])
        # read nc file
        nc=read_nc(dct['bc']['path'])
        
        # prepare df to store the gridded values of interest 
        nct=pd.DataFrame(columns=[tpl])
        # assign to the dataframe the tansposed matrix of
        # the corresponding value 
        nct[tpl]=nc.loc[tpl].transpose()
        # get list of areas selected by the user (defined at nuts3 level)
        arealist=list(pd.read_table(out_path+'nuts3_selection.txt', 
                          header=None, sep='\n')[0])    
        # macrosector name
        ms = dct_ms[tpl[1]]           
        for areait in arealist:
            df_areas = pd.DataFrame(area[nuts_lv]['area'].loc[areait])     
            df_areas['mult'] = df_areas['area'].multiply(nct.reindex(df_areas.index)[tpl])
            areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
            res['bc'].loc[areait]=areaxvar               
            res['delta'].loc[areait]=areaxvar*red[ms].loc[tpl[0]]/100  
        
        print('Saving results')
        res['per']=(res['delta'])/res['bc']*100
        res['value']=res['bc']-res['delta']
        res[['value', 'delta', 'per']].to_csv(out_path+nuts_lv, header=True, index=True, sep='\t', na_rep='NaN', mode='w')  

        res.index.rename('NUTS_Lv3', inplace=True)
        new=grd_int.reset_index().drop_duplicates(subset='NUTS_Lv3').set_index('NUTS_Lv3')
        new.drop(['ROW', 'COL', 'AREA_km2', 'POPULATION',  'CENTROID_X',  'CENTROID_Y'], axis=1, inplace=True)
        new['delta']=res['delta']
        
        nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']
        nuts_lvs.remove('NUTS_Lv3')
        
        for nuts_lv in nuts_lvs :
             # define df for results
            res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                               columns=['bc','delta'])
            nc=read_nc(dct['bc']['path'])
            # create a pd dataframe which has as columns the delta or base 
            # case value
            nct=pd.DataFrame(columns=[tpl])
            # assign to the dataframe the tansposed matrix of
            # the correscponding value 
            nct[tpl]=nc.loc[tpl].transpose()
            arealist = area[nuts_lv]['area'].index.levels[0]
            for areait in arealist:
                df_areas = pd.DataFrame(area[nuts_lv]['area'].loc[areait])     
                df_areas['mult'] = df_areas['area'].multiply(nct.reindex(df_areas.index)[tpl])
                areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
                res['bc'].loc[areait]=areaxvar               
            print('Saving results')
            res['delta']=new.groupby([nuts_lv])['delta'].sum()
            res['per']=(res['delta'])/res['bc']*100
            res['value']=res['bc']-res['delta']
            res[['value', 'delta', 'per']].to_csv(out_path+nuts_lv, header=True, index=True, sep='\t', na_rep='NaN', mode='w')  
                       

    else: 
        for nuts_lv in nuts_lvs:
            # prepare pd dataframe for results (index are the names of the 
            # geographical units in each level and keys are the delta and the base
            # case values)
            res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                               columns=['bc','delta'])
            # for the delta and base case value
            for key in ['bc','delta']:
                # read the corresonding nc
                nc = read_nc(dct[key]['path'])   
                # drop level if necessary (e.g. when reading delta_concentration)
                tpl = dct[key]['var']
                if len(tpl)!=nc.index.nlevels:
                       nc = nc.reset_index(level=0, drop=True)
    
                # create a pd dataframe which has as columns the delta or base 
                # case value
                nct=pd.DataFrame(columns=[tpl])
                # assign to the dataframe the tansposed matrix of
                # the correscponding value 
                nct[tpl]=nc.loc[tpl].transpose()
                arealist = area[nuts_lv]['area'].index.levels[0]
    
                aggr=ast.literal_eval(dct[key]['aggregation'])
                if len(aggr)==2:
                    opt1=aggr[0]
                    opt2=aggr[1]
                else:
                    opt1=aggr
                    opt2=None
                     
                for areait in arealist:
                # create a pd dataframe with the areas of each cell 
                # in the geographical unit
                    # if the aggregation mode is average (e.g. concentration levels)
                    if opt1 =='avg':
                        # area averaged values 
                        df_areas = pd.DataFrame(area[nuts_lv][opt2].loc[areait])
                        df_areas['mult'] = df_areas[opt2].multiply(nct.reindex(df_areas.index)[tpl])
                        areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()                                      
                        areatot=df_areas[df_areas['mult'].notnull()][opt2].sum()
                        # if the aggregation mode is sum (e.g. mortality)
                    elif opt1 =='sum':
                        # sum the values in each cell considering the fraction of the 
                        # cell belonging to a geographical entity
                        df_areas = pd.DataFrame(area[nuts_lv]['parea'].loc[areait])
                        df_areas['mult'] = df_areas['parea'].multiply(nct.reindex(df_areas.index)[tpl])
                        areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
                        areatot = 1
           
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
    
    pass  
#    aggrinp_txt='D:/programs/sherpa/app/data/temp/aggrinp.txt'
#    aggregationbyarea(aggrinp_txt)
#    
#    a= {
#        "delta": {
#            "path": "D:/programs/sherpa/app/data/temp/delta_concentration.nc",
#            "var": "delta_concentration",
#            "aggregation": "('avg','pop')"},
#    	"bc": {
#            "path": "D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_concentrations/BC_conc_PM25_Y.nc",
#            "var": "conc",
#            "aggregation": "('avg','pop')"},
#    	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
#    	"output-dir":"D:/programs/sherpa/app/data/temp/"
#        }
#    b= {
#        "delta": {
#            "path": "D:/programs/sherpa/app/data/temp/delta_concentration.nc",
#            "var": "delta_concentration",
#            "aggregation": "('avg','area')"},
#    	"bc": {
#            "path": "D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_concentrations/BC_conc_PM25_Y.nc",
#            "var": "conc",
#            "aggregation": "('avg','area')"},
#    	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
#    	"output-dir":"D:/programs/sherpa/app/data/temp/"
#        }
#        
#        
#    
#    c= {
#        "bc": {
#            "path": "D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_emissions/BC_emi_PM25_Y.nc",
#            "var": "('PPM','Nsnaps07')",
#            "aggregation": "'sume'"},
#    	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
#    	"output-dir":"D:/programs/sherpa/app/data/temp/"
#        }
#    
#    b= {
#    "delta": {
#        "path": "D:/programs/sherpa/app/data/temp/healthimp.nc",
#        "var": "d_mort",
#        "aggregation": "'sum'"},
#	"bc": {
#        "path": "D:/programs/sherpa/app/data/temp/healthimp.nc",
#        "var": "v_mort",
#        "aggregation": "'sum'"},
#	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
#	"output-dir":"D:/programs/sherpa/app/data/temp/"
#    }
#    