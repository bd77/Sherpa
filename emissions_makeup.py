import numpy as np 
import pandas as pd  
'''
NAME
    Import info on grid points attribution to nuts or specific area type from ascii file
PURPOSE 
    Import info on grid points attribution to nuts or specific area type from ascii file/s. 
    If the file is single then it must contain the column 'Area [km2]' relative to % of the area in the finest nut, 
    this datum will be set to each nut but it will then aggregated for larger nuts when nutsarea will be calculated
    If the files are two, then each nut will have its own % area for each grid point, then the data will be merged here
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''    
def read_nuts_area(filenuts,calcall=False,nullnut=None):
    nuts_def= filenuts +'.txt'
    nuts_info = pd.read_csv(nuts_def,delimiter="\t")
    nuts_info=nuts_info.dropna(axis=1,how='all')
    nutsnames=list(nuts_info.columns[~nuts_info.columns.isin(['POP','COL','ROW','Area [km2]','LAT','LON'])])
    #optional 'nut' comprising all grid points
    if calcall :
        #nutsnames.insert(0, 'ALL')
        nutsnames.insert(0, 'ALL_'+nutsnames[0])
        nuts_info[nutsnames[0]]=nutsnames[0] 
    nuts_info['grid']=['_'.join(str(i) for i in z) for z in zip(nuts_info['COL'],nuts_info['ROW'])]      
    if 'Area [km2]' in nuts_info.columns:
        nuts_area=pd.concat(map(lambda p: nuts_info['Area [km2]'],nutsnames),axis=1)
        #nuts_area.index=nuts_info['grid']
        nuts_area.columns=nutsnames  
       #nuts_info=nuts_info[nutsnames]
    else:
        sys.exit("missing infos on grid cells area per nut")

    #aggregate data for each nut, create a dictionary
    nuts_info_all={}
    for nut in nutsnames:
        #create a multindex
        index = pd.MultiIndex.from_tuples(list(zip(nuts_info[nut],nuts_info['grid'])), names=['nutname','grid'])
        nut_info=pd.Series(list(nuts_area[nut]), index=index)
        nut_info=nut_info.to_frame(name='area')
        #aggregate data on these nuts if necessary
        nut_info_nut=nut_info.groupby(level=[0,1]).sum()    
        #find total area
        grid_area_tot=nut_info_nut.groupby(level=['grid']).sum()
        nut_info_nut['parea']=nut_info_nut/grid_area_tot
        nut_info_nut.loc[nut_info_nut['area']==0,'parea']=0.
        #eventually remove the fillng code
        if nullnut is not None:
            nut_info_nut=nut_info_nut.drop(nullnut, level='nutname')
        nuts_info_all[nut]=nut_info_nut
            
    return nuts_info_all

'''
NAME
    string manitulation for ara names
PURPOSE
    remove strange characters and shorts the names if necessary
PROGRAMMER(S)
    Denise Pernigotti 
REVISION HISTORY
    
REFERENCES
    
'''
import re
def name_short(name,lmax=12):
    sep = '/'
    name = name.split(sep, 1)[0].strip()
    name=name.replace('BELGIQUE-BELGIE', "Belgium")
    #remove between parenthesis text
    name=re.sub("[\(\[].*?[\)\]]", "", name).strip()
    name=name.replace('UNITED', 'Un')
    names=name.title()
    pattern = re.compile(r'\W+')
    names=names.replace('Poranesna Jugoslovenska Republika Makedonija','J Makedonija')
    names=names.replace("Prov.", "Prv")
    names=names.replace("Region", "Reg")      
    names=names.replace('Republic', "Rep")
    name=name.replace('Republika', "Rep")
    names=names.replace('Kreisfreie Stadt','')
    names=names.replace('Deutschsprachige Gemeinschaft','')

    if len(names)>lmax:
        names = pattern.split(names)
        #remove short strings (coniugation)       
        names=' '.join([word for word in names if len(word) >2])
        names=names.replace('And', '')
        names=names.replace('North East', "NE")
        names=names.replace('North West', "NW")
        names=names.replace('South West', "SW")
        names=names.replace('Northern', "N")
        names=names.replace('North', "N")
        names=names.replace('South', "S")
        names=names.replace('East', "E")
        names=names.replace('West', "W")
        if len(names)>lmax:
            vowels = ('a', 'e', 'i', 'o', 'u')
            names=''.join([l for l in names if l not in vowels])
            if len(names)>lmax:
                names = pattern.split(names)  
                names=' '.join([word[:3] for word in names])
                if len(names)>lmax:
                    names=names[:lmax]
    names=names.strip()
    return names
'''
NAME
    Create the emissions reducted to a given area for a predefined amount for each precursor
PURPOSE
    Create the emissions reducted to a given area for a predefined amount (1=100%) for each precursor. 
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
def create_emissionsdelta(emissions,emissions_dprecursor,narea):
    emissionsdelta=emissions.multiply(emissions_dprecursor,level=0,axis=0)
    #create dataframe of zeroes
    emissionsdelta_area=pd.DataFrame(np.zeros(emissionsdelta.shape),columns=emissionsdelta.columns)
    #set only values in the area to non zero delta, multiplied for the relative area in the nut area
    emissionsdelta_area[narea.index]=emissionsdelta[narea.index].values*narea[narea.index].values
    emissionsdelta_area.index=emissionsdelta.index
    return emissionsdelta_area
'''
NAME
    Read information on areas of emissions and their metadata as prety names from ascii files
PURPOSE
    Read information on areas of emissions and their metadata as prety names from ascii files
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''    
def areas_setup(ascii_path,areaname,aggr_info):
    grid_txt=ascii_path + '\\grid_intersect'
    gcities_txt=ascii_path + '\\grid_int_gcities'
    fua_txt=ascii_path + '\\grid_int_fua'
    #list of true names for areas IDs
    codes_names=pd.Series({'NUTS_Lv0':'nuts0','NUTS_Lv1':'nuts1','NUTS_Lv2':'nuts2','NUTS_Lv3':'nuts3','FUA_CODE':'fua','GCITY_CODE':'gcities'})
    codes_txt={k: ascii_path + '\\' + codes_names[k] +'_names.txt'for k in codes_names.keys()}
    #list of points used as receptors. Only 'id','lon','lat' and 'Station name' are required columns 
    targets_txt=ascii_path + '\\' + areaname + '_targets.txt'

    #grab nuts/are info from txt
    nuts_info=read_nuts_area(grid_txt,calcall=True)
    nuts_info.update(read_nuts_area(gcities_txt,nullnut='LAND000'))
    nuts_info.update(read_nuts_area(fua_txt,nullnut='LAND000'))

    #grab true names for each area
    area_names_long={k: pd.read_csv(codes_txt[k],delimiter="\t",encoding = "ISO-8859-1",skiprows=1,index_col=0,header=None) for k in codes_names.keys()}
    #reduce string length if needed
    area_names={}
    for k in area_names_long.keys():
        area_names[k]=area_names_long[k][1].apply(name_short)
    area_names['ALL_NUTS_Lv0']=pd.Series({'ALL_NUTS_Lv0' : 'Europe', 'other' : 'other'},name='EU')
    #define the aggregation and increments calculation type depending on aggr_zones

    if aggr_info=='city':
        wantedorder=pd.Series( {'3' : 'GCITY_CODE', '2' : 'FUA_CODE', '1' :'NUTS_Lv0','0':'ALL_NUTS_Lv0'})
    elif aggr_info=='nuts':
        wantedorder=pd.Series( {'4' : 'NUTS_Lv3','3' : 'NUTS_Lv2', '2' : 'NUTS_Lv1', '1' :'NUTS_Lv0','0':'ALL_NUTS_Lv0'})
    outdic={'nuts_info':nuts_info, 'area_names':area_names,'wantedorder':wantedorder}
    return outdic

