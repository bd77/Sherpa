import numpy as np
import pandas as pd
'''
NAME
    Implementation of haversine formula (form lat lon to distances in km) for vectors
    Calculate the great-circle distance between two points on the Earth surface.
PURPOSE
    Implementation of haversine formula (form lat lon to distances in km) for vectors
    :input: one 2-tuples, and a vector of 2-tuples containing the latitude and longitude of a point
    in decimal degrees and a vector.
    Example: haversine((45.7597, 4.8422), (lat, lon))
    :output: Returns the distance in km bewteen the the point to all other points.PROGRAMMER(S)
    Denise Pernigotti from https://github.com/mapado/haversine/blob/master/haversine
REVISION HISTORY
    
REFERENCES
    
'''
AVG_EARTH_RADIUS = 6371  # in km


def haversine_vec(lon1,lat1,lon_vec2,lat_vec2):

    # calculate haversine
    dlat = np.radians(lat_vec2) - np.radians(lat1)
    dlon = np.radians(lon_vec2) - np.radians(lon1)
    d = np.sin(dlat * 0.5) ** 2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat_vec2)) * np.sin(dlon * 0.5) ** 2
    h = 2 * AVG_EARTH_RADIUS * np.arcsin(np.sqrt(d))
    return h # in kilometers
'''
NAME
    Given a grid point find whch areaid is pertaining (with the greater percentage)
PURPOSE
    Given a grid point find whch areaid is pertaining (with the greater percentage)
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
def find_target_info(areaid,areacells,x_y):
      area_find=areacells[areaid].swaplevel(i=-2, j=-1, axis=0)
      sorted_target_parea=area_find.loc[x_y].sort_values('parea',ascending=False)
      if len(sorted_target_parea)>0:
          target_area=list(sorted_target_parea.index)[0]
          #if sorted_target_parea['parea'].values[0] <0.3:
              #print ("for nut ",areaid," the tested point is for ","{:3.1f}".format(list(sorted_target_parea['parea'])[0]*100),"% in ",target_area, " calculation not performed")
          #else:
          if sorted_target_parea['parea'].values[0] >=0.3:
              #print ("for nut ",areaid," the tested point is for ","{:3.1f}".format(list(sorted_target_parea['parea'])[0]*100),"% in ",target_area)
            #select grid point in the target_area
              narea=areacells[areaid].loc[target_area]['parea']
            #a check to remoce spurious grid cells
              if len(narea[narea==0].index)>0:
                  #print("There are", len(narea[narea==0].index),"grid points with zero area in the selected area, they are removed")
                  narea=narea.loc[narea>0]
              area_cells=narea.index
            #print ('there are ',len(area_cells),' points with nut ',rads,'=',target_area)
              target_info=pd.Series([target_area,len(area_cells)],index=['areaid','ncells'],name=areaid)
              return target_info

'''
NAME
    Read user receptor list and check consistency
PURPOSE
    Read user receptor list and check consistency
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
''' 
def user_receptors(asciipath, testarea,coord,conc):     
    targets_txt=asciipath + '\\' + testarea + '_targets.txt'        
    #read list of receptors and check its consistency
    receptors=pd.read_csv(targets_txt,delimiter="\t",encoding = "ISO-8859-1")
    receptors.columns=[x.lower() for x in  receptors.columns]
    #check if information in receptors are ok and without duplicates
    required_info=['id','station name','lon','lat']
    required_info_present=list(filter(lambda x: x in receptors.columns, required_info))
    missing_info=list(set(required_info)-set(required_info_present))
    if len(required_info_present)<len(required_info):
        sys.exit("the receptors location file is missing some information on "+ missing_info)
    if len(receptors['id'].unique())<len(receptors['id']):
        sys.exit("there are duplicates in receptors ids")  
    if receptors['lon'].min()<coord['lon'].min() or receptors['lon'].max()>coord['lon'].max():
        sys.exit("some receptors have longitude outside the domain") 
    if receptors['lat'].min()<coord['lat'].min() or receptors['lat'].max()>coord['lat'].max():
        sys.exit("some receptors have latitude outside the domain")
    alldist_km=pd.concat(list(map(lambda st: haversine_vec(receptors.loc[st,'lon'],receptors.loc[st,'lat'],coord['lon'],coord['lat']),receptors.index)),axis=1)

    receptors['id']=receptors['id'].str.strip()
    receptors['target_idx']=alldist_km.idxmin()
    receptors['dist_km']=alldist_km.min()
    receptors['lon_grid']=pd.Series({st:coord.loc[receptors.loc[st,'target_idx'],'lon'] for st in receptors.index})
    receptors['lat_grid']=pd.Series({st:coord.loc[receptors.loc[st,'target_idx'],'lat'] for st in receptors.index})
    receptors['model_conc']=pd.Series(conc[receptors['target_idx']].values[0,])
    count_idx=receptors.pivot(columns='target_idx', values='target_idx').count()
    if count_idx.max()>1:
        print('There are duplicates in target_idx')
        print(count_idx.loc[count_idx>1,])
    a =count_idx.reset_index()
    a.rename(columns={0: 'duplicates'}, inplace=True)
    receptors=pd.merge(receptors,a)
    receptors.index=receptors['id']
    receptors.drop('id', axis=1, inplace=True)
    return receptors
