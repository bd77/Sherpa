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
import numpy as np
from time import time
import json
import ast 


#import os as os 

from sherpa_auxiliaries import read_nuts_area 
from sherpa_auxiliaries import read_nc

def module9_aggregation(aggrinp_txt):
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
        'NUTS0', 'NUTS1','NUTS2', 'NUTS3' 
        
    '''
#    aggrinp_txt='D:/programs/sherpa-v.2.0-beta.2/app/data/temp/aggregation2.txt'
    json_file = open(aggrinp_txt)
    json_str = json_file.read()
    dct = json.loads(json_str)
    grd_int_txt = dct['grid-intersect']
    out_path =  dct['output-dir']
    start = time()    # read the grid intersect (fua or nuts)    
    area=read_nuts_area(grd_int_txt, calcall=True)
    # levels (they are colled the same in fua and nuts)
    nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']

    dct_ms={'MS1': 'Nsnaps01',
     'MS10': 'Nsnaps10',
     'MS2': 'Nsnaps02',
     'MS3': 'Nsnaps03',
     'MS4': 'Nsnaps04',
     'MS5': 'Nsnaps05',
     'MS6': 'Nsnaps06',
     'MS7': 'Nsnaps07',
     'MS8': 'Nsnaps08',
     'MS9': 'Nsnaps09',
     'ALL': 'ALL'}
    inv_dct_ms = {v: k for k, v in dct_ms.items()}
    # @todo I have to try with U+00B5 and U+00B3
    dct_units={'conc': '[\u03bcg/m\u00B3]', '[delta_concentration]':'[\u03bcg/m\u00B3]',
               'v_dll_pp':"[dll/(person year)]",'d_dll_pp':'[dll/(person year)]', 
               'v_mort': '[people/year]', 'd_mort': '[people/year]'}
    grd_int=pd.read_csv(grd_int_txt+'.txt', sep='\t')
    # if we are aggregating emission reductions we need to consider the 
    # exact area defined by the user or reductions are smeared out on the grid. 
    if ast.literal_eval(dct['bc']['aggregation'])=='sume':
        nuts_lv = 'NUTS_Lv3'
        # get tuple defining precursor and sector
        tpl = ast.literal_eval(dct['bc']['var'])
        t=(tpl[0],dct_ms[tpl[1]])
        # read reduction file 
        red = pd.read_table(out_path+'user_reduction.txt', index_col=['POLL'])
        # define df for results

        
        # read nc file
        nc=read_nc(dct['bc']['path'])
        
        # prepare df to store the gridded values of interest 
       
        # get list of areas selected by the user (defined at nuts3 level)
        arealistall=list(pd.read_table(out_path+'nuts3_selection.txt', 
                          header=None, sep='\n')[0]) 
        # remove areas that are not in the domain 
        arealist = set(arealistall) - (set(arealistall) - set(grd_int[nuts_lv]))
        
        # macrosector name
        ind_areas = area[nuts_lv]['area'].index.levels[0]
        if t[1]!='ALL':
            ms_list=[t[1]]

        elif t[1]=='ALL': 
            mss=list(dct_ms.values())
            mss.remove('ALL')
            ms_list=mss
            
        for ms in ms_list:
            tpl_new=(t[0], ms)
            nct=pd.DataFrame(columns=[(tpl_new)])
            # assign to the dataframe the tansposed matrix of
            # the corresponding value 
            nct[tpl_new]=nc.loc[tpl_new].transpose()
            ind_ms=[ms]*len(ind_areas)
            ind=pd.MultiIndex.from_tuples(list(zip(ind_ms,ind_areas)), names=('sector', 'area'))
            res = pd.DataFrame(index=ind,
                           columns=['bc','delta'])
            for areait in arealist:
                df_areas = pd.DataFrame(area[nuts_lv]['area'].loc[areait])     
                df_areas['mult'] = df_areas['area'].multiply(nct.reindex(df_areas.index)[tpl_new])
                areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
                res['bc'].loc[ms,areait]=areaxvar               
                res['delta'].loc[ms,areait]=areaxvar*red[inv_dct_ms[tpl_new[1]]].loc[tpl_new[0]]/100  
        
        print('Saving results')
        if t[1]!='ALL':
            res=res.loc[t[1]]
        else:
            res=res.sum(level=[1])
            
        res['per']=(res['delta'])/res['bc']*100
        res['value']=res['bc']-res['delta']
       
        res[['value', 'delta','per']].rename(
                columns={'value':'value[Mg]', 'delta':'delta[Mg]', 'per':'per[%]'}).to_csv(
                        out_path+nuts_lv[0:4]+nuts_lv[-1]+'.txt', header=True, index=True, sep='\t', na_rep='NaN', mode='w',encoding='utf-8')  

        res.index.rename('NUTS_Lv3', inplace=True)
        new=grd_int.reset_index().drop_duplicates(subset='NUTS_Lv3').set_index('NUTS_Lv3')
        new.drop(['ROW', 'COL', 'AREA_km2', 'POPULATION',  'CENTROID_X',  'CENTROID_Y'], axis=1, inplace=True)
        new['delta']=res['delta']
        new['bc']=res['bc']
        
#        nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']
        nuts_lvs.remove('NUTS_Lv3')
        
        for nuts_lv in nuts_lvs :

             # define df for results
            res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                               columns=['bc','delta'])
             
            print('Saving results')
            res['bc']=new.groupby([nuts_lv])['bc'].sum()
            res['delta']=new.groupby([nuts_lv])['delta'].sum()
            res['per']=(res['delta'])/res['bc']*100
            res['value']=res['bc']-res['delta']
            res[['value', 'delta','per']].rename(
                columns={'value':'value[Mg]', 'delta':'delta[Mg]', 'per':'per[%]'}).to_csv(
                        out_path+nuts_lv[0:4]+nuts_lv[-1]+'.txt', header=True, index=True, sep='\t', na_rep='NaN', mode='w',encoding='utf-8')  
                       

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
                    # shape file which are outside the domain 
                    if areatot is not 0: 
                        value=areaxvar/areatot
                    else:
                        value= float('NaN')
                    res[key].loc[areait]=value           
            
            units = dct_units[dct['bc']['var']]
            print('Saving results')
            res['per']=(res['delta'])/res['bc']*100
            res['value']=res['bc']-res['delta']
            res[['value', 'delta', 'per']].rename(
                columns={'value':'value'+units, 'delta':'delta'+units, 'per':'per[%]'}).to_csv(
                        out_path+nuts_lv[0:4]+nuts_lv[-1]+'.txt', header=True, index=True, sep='\t', na_rep='NaN', mode='w', 
                        encoding='utf-8')  
        end = time()
        print('Calculation time  for aggregation', end-start)
    
if __name__ == '__main__': 
    
    pass  
