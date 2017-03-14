import numpy as np
import pandas as pd
'''
NAME
    aggregate snaps as requred by e-rep 
PURPOSE
    aggregate snaps as requred by e-rep
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''  
def dc_snapaggregate(alldc_snaps):
    #aggregate sources as required in e-rep
    sources=pd.Series(['Traffic','Industry','Industry','Industry','Agriculture','Residential','Offroad','Other','Other','Other'],
                        index=['Nsnaps7','Nsnaps1','Nsnaps3','Nsnaps4','Nsnaps10','Nsnaps2','Nsnaps8','Nsnaps5','Nsnaps6','Nsnaps9'])
    alldc=alldc_snaps
    alldc['source']=sources.loc[alldc_snaps.index]
    alldc.set_index('source', append=True, inplace=True) 
    alldc=alldc_snaps.groupby(level=1).sum()
    return alldc

'''
NAME
    calculate increments 
PURPOSE
    calculate increments
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''  
def dc_increments(alldc,aggr_zones='city'):    #calculate increments
    alldc_inc={}
    alldc_inc['ALL_NUTS_Lv0']=alldc['ALL_NUTS_Lv0']-alldc['NUTS_Lv0'] 
    if aggr_zones=='city': 
        if 'FUA_CODE' in alldc.columns:
            alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']-alldc['FUA_CODE']
            if 'GCITY_CODE' in alldc.columns: 
                alldc_inc['FUA_CODE']=alldc['FUA_CODE']-alldc['GCITY_CODE']
                alldc_inc['GCITY_CODE']=alldc['GCITY_CODE']
            else:
                alldc_inc['FUA_CODE']=alldc['FUA_CODE']               
        else:
           alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']
    elif aggr_zones=='nuts':
        alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']-alldc['NUTS_Lv1']
        if 'NUTS_Lv2' in alldc.columns:
            alldc_inc['NUTS_Lv1']=alldc['NUTS_Lv1']-alldc['NUTS_Lv2']
            if 'NUTS_Lv3' in alldc.columns: 
                alldc_inc['NUTS_Lv2']=alldc['NUTS_Lv2']-alldc['NUTS_Lv3']
                alldc_inc['NUTS_Lv3']=alldc['NUTS_Lv3']
            else:
                alldc_inc['NUTS_Lv2']=alldc['NUTS_Lv2']               
        else:
           alldc_inc['NUTS_Lv1']=alldc['NUTS_Lv1']
     
    alldc_inc=pd.DataFrame.from_dict(alldc_inc)
    return alldc_inc        
'''
NAME
    Summation of dc per area and eventually per precursor
PURPOSE
    Summation of dc per area and eventually per precursor
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''  
def dc_areasum(dc_all,narea,liv=1):
    dc_tot=pd.concat(list(map(lambda areaname:(dc_all[narea[areaname].index]*narea[areaname].values).sum(axis=1) ,narea.keys())),axis=1)
    dc_tot.columns=list(narea.keys())
    #overall sum including sum over snaps
    if isinstance(dc_all.index, pd.core.index.MultiIndex):
        alldc_snaps=dc_tot.groupby(level=[liv]).sum()
    else:
        alldc_snaps=dc_tot
    return alldc_snaps
