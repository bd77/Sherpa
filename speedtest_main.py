    # main 
if __name__ == '__main__':
    
    import time
    import sys
    import os
    import numpy as np
    import pandas as pd  
     #import openpyxl
    from scipy.spatial import distance
    
    import sherpadf as sh
    
    #directory for emissions and SR relationship only
    path_nc='C:\\Users\\pernide\\AppData\\Local\\Sherpa\\app\\data\\input\\'
    path_txt='E:\\BKUP\\JRC\\E-REP\\input'
    ############################################### user input data
    pollname='PM25'
    test_all=False
     
    #load area information
    nuts_info=sh.read_nuts_area(path_txt + '\\selection\\grid_intersect')
    
    #load input files
    land_emissions_points=nuts_info['NUTS_Lv0'].index.get_level_values(1)
    ncinput=sh.clean_input(path_nc,pollname,emitting_points=land_emissions_points,include_natural=False)   
    for k, v in ncinput.items():
        locals()[k]=v
    
    test_number=20
    list_receptors=model.columns[np.arange(0,test_number)]
    emissions_idx=emissions_array['dimnames']['emissions_idx']
    snaps=emissions_array['dimnames']['snaps']
    precursors=emissions_array['dimnames']['precursors']
    start =time.perf_counter()    
    distall = distance.cdist(coordinates.loc[list_receptors,['x','y']], coordinates.loc[emissions_idx,['x','y']],metric='euclidean')
    end = time.perf_counter()
    print ("time elapsed dists on all points",end-start)
    for ix,idx in enumerate(list_receptors):
        dists_idx=pd.Series(distall[ix,:])
        dists_idx.index=emissions_idx
 
        start =time.perf_counter()  
        dc_orig=pd.concat(list(map(lambda p: sh.sherpa_df(p, model[idx],dists_idx,emissions,inner_radius),precursors))) 
        dc_orig.index=emissions.index
        end = time.perf_counter()
        print ("time elapsed orig erep with flatweight",end-start)
        
        start =time.perf_counter()    
        #Core of the sharpa calculations: calculate dc due to each gridpoint for each precursor and macrosector
        dc_orig_test=pd.concat(list(map(lambda p: sh.sherpa_df(p, model[idx],dists_idx,emissions),precursors))) 
        dc_orig_test.index=emissions.index
        end = time.perf_counter()
        print ("time elapsed orig erep without flatweight",end-start)
        dc_orig_array=np.zeros(shape=(len(snaps),len(emissions_idx),len(precursors)))
        for isp,prec in enumerate(precursors):
            dc_orig_array[:,:,isp]=np.array(dc_orig_test.loc[(prec,slice(None)),:])  

        start =time.perf_counter()
        dc_vec=sh.sherpa_np(model[idx],dists_idx,emissions_array) 
        end = time.perf_counter()
        print ("time elapsed array",end-start)
       
        start =time.perf_counter()
        dc_ne=sh.sherpa_ne(model[idx],distall[ix,:],emissions_array) 
        end = time.perf_counter()
        print ("time elapsed ne",end-start)
        
        if np.isclose(dc_vec,dc_ne).sum()<len(snaps)*len(precursors)*len(emissions_idx):
            print('******WARNING******** dc_ne and dc_vec are NOT close')
        if np.isclose(dc_orig_array,dc_ne).sum()<len(snaps)*len(precursors)*len(emissions_idx):
            print('******WARNING********  dc_ne and dc_orig_test are NOT close')
        

    if test_all:    
        dc_sum={}
        idx_list={}
        max_receptors=10000 #under the memory limit of the machine
        nprocs=int(len(model.columns)/max_receptors)+1
        for ip in np.arange(0,nprocs):
            minlist=ip*max_receptors
            maxlist=(ip+1)*max_receptors
            if maxlist > len(model.columns):maxlist=-1 
            idx_list['p'+str(ip)]=model.columns[minlist:maxlist]
        startall =time.perf_counter()
        for process in idx_list.keys():
            start =time.perf_counter()  
            distall = distance.cdist(coordinates.loc[idx_list[process],['x','y']], coordinates.loc[emissions_idx,['x','y']],metric='euclidean')
            end_dist=time.perf_counter()
            print ("time elapsed dist",end_dist-start)
            for ix,idx in enumerate(idx_list[process]):
                dc_ne=sh.sherpa_ne(model[idx],distall[ix,:],emissions_array) 
                dc_sum[idx]=dc_ne.sum(axis=(1,2))
            end = time.perf_counter()
            print ("time elapsed process",end-start)
    
        endall = time.perf_counter()
        print ("time elapsed ne on all points",endall-startall)
        dc_sum_df=pd.DataFrame.from_dict(dc_sum,orient='index')
        dc_sum_df.to_csv(path_or_buf='dc_sum')
    
           
    
    
