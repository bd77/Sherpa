# main 
if __name__ == '__main__':
    
    import time
    import sys
    import os
    import numpy as np
    import pandas as pd  
    from mpl_toolkits.basemap import Basemap    
    import matplotlib.pyplot as plt
    from scipy.spatial import distance
    
    import sherpadf as sh
    
    #directory for emissions and SR relationship only
    path_nc='C:\\Users\\pernide\\AppData\\Local\\Sherpa\\app\\data\\input\\'
    path_txt='E:\\BKUP\\JRC\\E-REP\\input'
    path_out='E:\\BKUP\\JRC\\E-REP\\output'
    ############################################### user input data
    pollname='PM25' #may be '25' or '10'
    testarea='Paris' #may be any area as long as the file testarea_targets.txt is present in input, contains a list of lat/lon
    aggr_zones='nuts' #may be 'city' or 'nuts' 
    include_natural=False #include or not natiral PM concentration
    exclude_sea_emissions=True #Emissions on sea are not considered in the present version of sherpa
    outfig='png' #'pnd' of 'pdf'
    ############################################### 
    ############################################### input files

 
    #If not present create the output/testarea directory
    if not os.path.exists(path_out +'\\'+testarea):
        os.makedirs(path_out +'\\'+testarea)
    
    #load area information
    areainfo=sh.areas_setup(path_txt+'\\selection',testarea,aggr_zones)
    for k, v in areainfo.items():
        locals()[k]=v
    
    #load input files
    land_emissions_points=[]
    if exclude_sea_emissions:
        land_emissions_points=nuts_info['NUTS_Lv0'].index.get_level_values(1)
    ncinput=sh.clean_input(path_nc,pollname,emitting_points=land_emissions_points,include_natural=False)   
    for k, v in ncinput.items():
        locals()[k]=v
 
    receptors=sh.user_receptors(path_txt,testarea,coordinates,concentration)
    count_idx=receptors.pivot(columns='target_idx', values='target_idx').count()  
    precursors=model.index.get_level_values(1).unique()
    #emissions=create_emissionsdelta(emissions,emissions_dprecursor,narea_all) 
     
    #test_point={'lon':8.8125,'lat':45.84375}
    #test_point={'lon':9.191383,'lat':45.464211} #Milan
    #test_point={'lon':9.1875,'lat':45.46875} #closer grid point of Milan
    #test_point={'lon':2.349014,'lat':48.864716} #Paris 
    #find distance from all other points in emissions data in km
    #emissions['dist_greatc']=mapnu(lambda x:great_circle((test_point['lat'],test_point['lon']),(emissions.loc[x,'lat'],emissions.loc[x,'lon'])).km,emissions.index)
 
        
    dc_inc_all={}
    dc={}
    dc_ppm={}
    target_allinfo={}
    for idx in count_idx.index:
        #For the selected point find all information on area  ids
        target_info=pd.concat(list(map(lambda areaid: sh.find_target_info(areaid,nuts_info,idx),nuts_info.keys())),axis=1).transpose()
        target_info['areaname']=pd.Series(dict(zip(target_info.index,map(lambda x: ''.join(list(area_names[x].loc[target_info.loc[x,'areaid']])).title(),target_info.index))))     
        #select grid points in the same area as the target, with their percentual area    
        narea=pd.concat(list(map(lambda areaname:nuts_info[areaname].loc[target_info.loc[areaname,'areaid']]['parea'],target_info.index)),axis=1)
        narea.columns=target_info.index
        narea=narea.fillna(0)
        #narea = dict(zip(target_info.index, narea))          
        #create an emissions delta with zeroes out of the area and emissions detuced by emissions_delta multiplies for the realtive grid area inside the area  

        #Core of the sharpa calculations: calculate dc due to each gridpoint for each precursor and macrosector
        start =time.perf_counter()    
        dists_idx=pd.Series(np.ndarray.flatten(distance.cdist(coordinates.loc[idx,['x','y']][np.newaxis, :],coordinates[['x','y']], metric='euclidean')),index=coordinates.index)
        dc[idx]=pd.concat(list(map(lambda p: sh.sherpa_df(p, model[idx],dists_idx,emissions,inner_radius),precursors))) 
        dc[idx].index=emissions.index
        end = time.perf_counter()
        #print ("time elapsed ",end-start)

        #aggregate dc per precursor, area, calculate increments and relative values
        alldc=sh.dc_snapaggregate(sh.dc_areasum(dc[idx],narea))
        dc_inc=sh.dc_increments(alldc,aggr_zones)*100./concentration[idx].values[0] 
        area_present=dc_inc.columns
        dc_inc['total']= dc_inc.sum(axis=1)
        
        #aggregate results per precursor
        alldc_prec=sh.dc_areasum(dc[idx],narea,liv=0)
        dc_inc_p=sh.dc_increments(alldc_prec,aggr_zones)*100./concentration[idx].values[0]
        dc_inc_p['total']= dc_inc_p.sum(axis=1)
   
        wantedorder_present=pd.Series(list(filter(lambda x: x in area_present, wantedorder))).to_frame(name='areaid')
        wantedorder_present['areaname']=pd.Series(list(target_info.loc[wantedorder_present['areaid'],'areaname']))
        #check for duplicated names in nuts
        if 'NUTS_Lv1' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv1','areaname']='1_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv1','areaname']
        if 'NUTS_Lv2' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv2','areaname']='2_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv2','areaname']

         #rename fue elements (same names as city core)
        if 'FUA_CODE' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']='fua_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']

        #build a dataframe with constant values in the smallest areas (excluding 'ALL_NUTS_Lv0' and 'NUTS_Lv0')
        smallareas=list(set(area_present)-set(['ALL_NUTS_Lv0','NUTS_Lv0']))
        #avoid double counting of grid increments
        narea_inc=sh.dc_increments(narea,aggr_zones)
        
        if len(smallareas)>0:
            dc_inc_flat= pd.concat(list(map(lambda p: narea_inc[p]*(dc_inc[p].sum()),smallareas)),axis=1).sum(axis=1) 
        dc_inc_flat= dc_inc_flat.reindex(index=dc[idx].columns) 
        
        #set up appropriate names for the areas
        wantedorder_present=wantedorder_present.append(pd.Series({'areaid':'total', 'areaname':'total'}), ignore_index=True)
        dc_inc=dc_inc[wantedorder_present['areaid']]
        #add natural sources
        if include_natural:
            pnsize=pollname[2:4]
            natural=pd.DataFrame(0, index=['Salt','Dust'], columns=dc_inc.columns)
            natural.loc['Dust','total']=dust.loc['pDUST-'+pmsize,idx].values*100./concentration[idx].values[0]
            natural.loc['Salt','total']=salt.loc['pSALT-'+pmsize,idx].values*100./concentration[idx].values[0]
            dc_inc=dc_inc.append(natural)
        #plots
        fig={}
        fig[1]=sh.plot_dict(dc_inc_flat,idx,coordinates['lat'],coordinates['lon'])
        plt.close('all')
        fig[2]=sh.plot_bar(dc_inc,wantedorder_present)
        plt.close('all')
        fig[3]=sh.plot_bar(dc_inc_p,wantedorder_present)
        plt.close('all')
        dc_inc.columns=wantedorder_present['areaname']
  
        for ids in list(receptors[receptors['target_idx']==idx].index):
          fig[1].savefig(path_out +'\\'+testarea+'\\'+ids+'_'+pollname+'_'+aggr_zones+'_sec_map.'+outfig)
          fig[2].savefig(path_out +'\\'+testarea+'\\'+ids+'_'+pollname+'_'+aggr_zones+'_sec_bars.'+outfig)
          fig[3].savefig(path_out +'\\'+testarea+'\\'+ids+'_'+pollname+'_'+aggr_zones+'_prec_bars.'+outfig)
          dc_inc.to_html(path_out +'\\'+testarea+'\\'+ids+'_'+pollname+'_'+aggr_zones+'_total_table.html',classes='table')
          dc_inc_all[ids]=dc_inc.transpose()
          target_allinfo[ids]=target_info.transpose()
    #summarize info on grid points
    reform = {(outerKey, innerKey): values for outerKey, innerDict in target_allinfo.items() for innerKey, values in innerDict.items()}
    target_allinfo=pd.DataFrame(reform).transpose() 
    a =target_allinfo.reset_index()
    a.rename(columns={'level_0': 'id', 'level_1': 'area'}, inplace=True)
    b =receptors[['station name','target_idx','duplicates','model_conc','lon','lon_grid','lat','lat_grid','dist_km']].reset_index()
    summary= pd.merge(b, a)
    summary = summary.set_index(['id','area'])
    summary.to_html(path_out +'\\'+testarea+'\\AAA_summary_info.html',classes='table')

    reform = {(outerKey, innerKey): values for outerKey, innerDict in dc_inc_all.items() for innerKey, values in innerDict.items()}
    dc_inc_all=pd.DataFrame(reform).transpose()  
    #check total totals of explained mass on each receptor
    total_pm=pd.concat([receptors[['station name','lon','lat','target_idx']],dc_inc_all.groupby(level=[0]).sum().loc[receptors.index,]],axis=1)
    #COUNTaverage source contribution and its variability among receptors
    summary_src={}
    summary_src['count']=dc_inc_all.groupby(level=[1]).count().transpose()
    summary_src['mean']=dc_inc_all.groupby(level=[1]).mean().transpose()
    summary_src['sd']=dc_inc_all.groupby(level=[1]).std().transpose()
    reform = {(innerKey,outerKey): values for outerKey, innerDict in summary_src.items() for innerKey, values in innerDict.items()}
    summary_src=pd.DataFrame(reform).transpose()  
     # pdf from htl template as in http://stackoverflow.com/questions/27387923/combine-existing-figures-into-one-pdf-of-figure-python
    #path_wkthmltopdf = r'C:\Program Files\wkhtmltopdf\bin\wkhtmltopdf.exe'
    #config = pdfkit.configuration(wkhtmltopdf=path_wkthmltopdf)
    #pdfkit.from_file('output/template.html', 'output/output.pdf',configuration=config)

    
    
    
