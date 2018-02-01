'''
Created on Jun 23, 2015

Module description. This is the adaptation of the original module7 
developed by Denise Pernigotti. The main differences are the format of the 
input file (the grid_intersect, which had to be updated to TerraAria 
requirements) and the output figures which are now equivalent to the ones 
generated in the Atlas on Urban PM2.5. 

WARNING: I preferred keeping as much as the original code of Denise as 
possible, because it may be useful in the future. 

Inputs: 
    - emissions_nc: path to emissions nc file
    - concentration_nc: path to concentration nc file
    - natural_dir: path to the directory of the natural concentrations: 
        e.g. 
        natural_dir: './input/pDUST-pSALT/'
        which contains 4 files named as follows: pDUST-25-basecase.nc etc.
    - model_nc: path to model nc file
    - fua_intersect_dir: path to the directory where the grid_intersect.txt of fuas is located, 
                         e.g.:
                         fua_intersect_dir = './input/selection/gridnew/fua/'
                         fua_intersect_dir = './input/models/chimere_7km_2016_v1.7_fua/selection/'
    - dbf_dir: path to the directory where the shape of fua and nuts are located, 
                         for example:
                         e.g.
                         dbf_dir: 'D:/sherpa.git/Sherpa/input/selection/gridnew/'
                         dbf_dir: 'C:/Users/peduzem/AppData/Local/Sherpa/app/data/shapes'
    - targets_txt: path to the txt file where the targets coodrdinates are stored. 
    min example text file: 
    id	Station name	lon	lat
    UK09999	Liverpool	-2.9375	53.40625
    - outdir: directory to save results
    - aggr_zones: fua, nuts and fuaonly for the moment 
    - pollutant: 'PM25' or 'PM10' or 'NOx'
    - outfig: extension of figures to save, 'png' 
    - normalize=True Normalize figure values to 100%

output: - polar polots for emissions 
        - bar plot for concentrations 
        - legend plot (combining both legends)
        - excel files with summary of results

@author: Denise Pernigotti and Emanuela Peduzzi (EPE)
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

from matplotlib.ticker import AutoMinorLocator 
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from scipy.spatial import distance
from simpledbf import Dbf5


def figsize(scale):
    '''Squared figure 
    @author:EPE
    '''
    fig_width_pt = 450  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
    fig_height = fig_width #* golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size

def figsizer(scale): 
    '''Rectangualar figure 
    @author:EPE
    '''
    fig_width_pt = 390  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
    fig_height = fig_width * golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size



def read_nc(nc_file):
    '''
    NAME
        Reads SHERPA ncdf file with Python
    PURPOSE
        To read matrix data and put them in a multindexed dataframe
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
    
    REFERENCES
    
    '''
    nc_data = Dataset(nc_file, 'r')
    nc_dims = [dim for dim in nc_data.dimensions]
    nc_vars = [var for var in nc_data.variables]
    
    #sometimes the latitude is written just with lat as in model data
    latname=list(filter(lambda x: x in nc_vars, ['Lat','lat','latitude']))[0]
    lats = nc_data.variables[latname][:]
    lonname=list(filter(lambda x: x in nc_vars, ['Lon','lon','longitude']))[0]
    lons = nc_data.variables[lonname][:]
    
    #if there are three dimensional arrays
    if len(nc_dims)==3:
        ncz=str(list(set(nc_dims)-set(['latitude','longitude']))[0])
        nz=range(len(nc_data.dimensions[ncz]))
        if ncz=='pollutant':
            strpoll=nc_data.Order_Pollutant
            nznames=strpoll.split(', ')
        else:
            nznames=[ncz +"{:02d}".format(x+1) for x in nz]

    #create an index with lat and lon
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

    reform = {(outerKey, innerKey): values for outerKey, innerDict in allvar.items() for innerKey, values in innerDict.items()}
    df=pd.DataFrame(reform)
    return df.transpose()

def name_short(name,lmax=12):
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

def write_dict_nc(dc_dic,lats,lons,unit,filenc):
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

def plot_bar(dfdata, varplot, totalname, plot_opt='perc',
             x_label='Percentage of total mass', leg_loc='lower right', ftx=10, normalize=True):
    '''
    NAME
       Stacked barplot using plt
    PURPOSE
        Stacked barplot using plt,
        http://chrisalbon.com/python/matplotlib_stacked_bar_plot.html
        Input parameter ftx is the font size of labels, 16 is the default
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
        Revised by @author:EPE to make it consistent with the atlas
    REFERENCES
    '''   
    
    dfdata = dfdata[varplot['areaid']].transpose()
    addsum = dfdata.sum(axis=1)
    
    # Drop levels with zero values (we have to do this for cities like 
    # Liverpool that do not have commuting zones). I also remove levels whose 
    # emissions are really small and would not display well in the figures
    for key in addsum.index:
        if addsum[key]<=1e-4:
            print('Removing level ', key, ' because sum of contributions is small',addsum[key], '<1e-4' )
            dfdata.drop([key], inplace=True)
            varplot = varplot[varplot.areaid!=key]

    varplot_left = pd.Series(np.zeros(len(dfdata.index)), index=dfdata.index)    
    if addsum.loc['Total'] > 100 and normalize==True:
        print('WARNING: rescaling to 100')
        dfdata_plot = dfdata * 100 /  addsum.loc['Total']
    else:
        ncontrol_fill = 100 - addsum.loc['Total']
        dfdata_plot = dfdata
        dfdata_plot['External'] = 0
        dfdata_plot.loc['Total', 'External'] = ncontrol_fill
                   
    varnames = dfdata_plot.columns
    for ind in range(1, len(dfdata_plot.index)):
        varplot_left[ind] = varplot_left[ind-1]+addsum[dfdata_plot.index[ind-1]]
    if totalname in varplot['areaid'].values:
        varplot_left[totalname] = 0.

    #Colors ATLAS
    colors = {'Transport': 'red', 'Energy': 'blue',
          'Industry': 'gold', 'Production': '#8B4789', 'Waste': '#00FFFF',
          'Agriculture': 'green', 'Residential': '#0080FF', #'Agriculture': '#33FF99'
          'Offroad': '#8470ff', 'Extraction': '#00FF00',
          'Other': '#ee82ee', 'Natural': '#ffe7ba', 'Salt': '#ccffe5', 'Dust': '#ffffcc',
          'External': '#cdcdb4', 'bottom': 'None'} 
  
    yaxisnames = {'CCITY_CODE': 'City',
                  'FUA_CODE': 'Commuting Zone',
                  'NUTS_Lv1': 'NUTS_Lv1', 
                  'NUTS_Lv2': 'NUTS_Lv2',
                  'NUTS_Lv3': 'NUTS_Lv3',
                  'NUTS_Lv0': 'Rest of the country',
                  'ALL_NUTS_Lv0': 'Transboundary', 
                  'Total': 'Total'}

    # EPE: Display only greater city when the city core si below 300 km2
    # and aggregation is done only on FUAs
    yser = [yaxisnames[k] for k in varplot['areaid'].values]
    varplot = varplot.copy()
    varplot.loc[:,'yname']=yser
    areanames = varplot['yname']
    areanames.index = varplot['areaid']     
    if 'FUA_CODE' in areanames.index.values:
        if 'CCITY_CODE' not in areanames.index.values:
            areanames.loc['FUA_CODE'] = 'Greater city' 
    
    # Create the general blog and the "subplots" i.e. the bars
    plt.close('all')
    f, ax1 = plt.subplots(1, figsize=figsizer(0.9))
    if plot_opt == 'perc':
        ax1.set_xlim([0, 100])
    else:
        ax1.set_xlim([0, addsum.sum()])
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(ftx)
    # Set the bar width
    bar_width = 0.75

    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(varplot))]

    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2)-0.3 for i in bar_l]  # modified by EPI

    # Create a bar plot, in position bar_1
    ax1.barh(bar_l,
             # using the pre_score data
             dfdata_plot[varnames[0]],
             # set the width
             height=bar_width,
             # with pre_score and mid_score on the bottom
             left=list(varplot_left),
             label=varnames[0],
             # with color
             color=colors[varnames[0]])

    # Create a bar plot, in position bar_1
    for ivar in range(1, len(varnames)):
        ax1.barh(bar_l,
                 # using the mid_score data
                 dfdata_plot[varnames[ivar]],
                 # set the width
                 height=bar_width,
                 # with pre_score and mid_score on the bottom
                 left=list(varplot_left+dfdata_plot[
                           varnames[range(0, ivar)]].sum(axis=1)),
                 # with the label post score
                 label=varnames[ivar],
                 # with color
                 color=colors[varnames[ivar]])

    # set the y ticks with names
    plt.yticks(tick_pos, areanames, size=ftx, rotation=0)

    # set minor x ticks - added EMA
    minorLocator = AutoMinorLocator()
    ax1.xaxis.set_minor_locator(minorLocator)

    # Set the label and legends
    ax1.set_xlabel(x_label, fontsize=ftx)

    # Set a buffer around the edge
    plt.ylim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
    
    # Place initial letter on top bar    
    ypos = len(dfdata_plot.index.values)
    xposcum = 0          
    for letlabind in (np.arange(len(dfdata_plot.columns))):
        xpos = xposcum + dfdata_plot.loc['Total'][letlabind]/2
        xposcum = xposcum + dfdata_plot.loc['Total'][letlabind]
        letter = dfdata_plot.columns[letlabind][0]
        # place latter only if there is enough space
        if dfdata_plot.loc['Total'][letlabind] >= 5 and np.isnan(xpos) == False:
            ax1.text(xpos, ypos, letter, va= 'center',  ha= 'center', fontsize=(ftx+2))
    plt.close()
#    plt.show()
    return f

def plot_polar(emi_sum, prec, wantedorder_present, ftx=8):  
    """Polar plot for emissions, only if the point is in the city
    @author: peduzem

    """  
    titles_dct = {'NH3': "NH$_\mathbf{3}$", 
                  'SOx': "SO$_\mathbf{2}$", 
                  'NMVOC': "NMVOC", 
                  'PPM': 'PPM$_\mathbf{2.5}$', 
                  'NOx': "NO$_\mathbf{x}$"}
    
    dct_colors = {'CCITY_CODE': ['red', 'red', 2*'///', 'City'],
                      'FUA_CODE': ['blue', 'None', '', 'Greater city']}  
    
    sect_display = list(set(emi_sum.index.get_level_values(1)))
    sect_display.sort() 
    index = [sect_display[2], sect_display[4], sect_display[0], sect_display[1], sect_display[3]]

    # Polar plot are created only for cities (receptor point at least in FUA): 
    # initialize array for results: 
    df = []  
    # initialize array for max value 
    mv = []
    for key in dct_colors.keys():
        if key in wantedorder_present['areaid'].values:
            mv.append(max(emi_sum[key].loc[prec].reindex(index)/1000))
            df.append([emi_sum[key].loc[prec].reindex(index)/1000])
                  
    fig = plt.figure()
    # create figure
    fig = plt.figure(figsize=figsize(0.3), dpi=1000)
    ax = fig.add_subplot(111, projection="polar")
    ax.grid(True)
    ax.yaxis.grid(color='#aab0b7', lw=0.7)
    ax.xaxis.grid(color='#aab0b7', lw=0.7)
    # set border of figure area
    for spine in ax.spines.values():
        spine.set_edgecolor('#aab0b7')
        spine.set_zorder(0)
        spine.set_linewidth(1)
            
    # Define angles of y axis in polar plots
    theta = np.arange(len(index))/float(len(index))*2.*np.pi
    ax.set_xticks(theta)
    
    # Define tick labels for y axis
    for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_horizontalalignment('center')
            tick.label1.set_verticalalignment('top')
    # get radial labels away from plotted line (set angle)
    ax.set_rlabel_position(90)
        
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_verticalalignment('top')
    
    # make plots
    plots = []
    for it, dfdata in enumerate(df):
        line, = ax.plot(theta, df[it][0], color=dct_colors[df[it][0].name][0],
                        marker=None, label=None, zorder=3)
#        ax.fill(theta, df[it][0], color=dct_colors[df[it][0].name][1], alpha=0.3,
#                zorder=3, lw=0.3)
        ax.fill(theta, df[it][0], edgecolor=dct_colors[df[it][0].name][1],
                alpha=1, hatch=dct_colors[df[it][0].name][2],
                fill=False, lw=0.4, zorder=3)
        plots.append(line,)


    def _closeline(line):
        x, y = line.get_data()
        x = np.concatenate((x, [x[0]]))
        y = np.concatenate((y, [y[0]]))
        line.set_data(x, y)
    [_closeline(l) for l in plots]
    
    maxy = max(mv)
    # set yticks so that I have 5 y ticks for each figure
    ylim_lst = [1, 2, 5, 10, 15, 20, 40, 50, 60, 80, 100, 150, 200, 300, 400, 1000]
    for lim in ylim_lst:
        if maxy <= lim:
            ylim = lim
            ax.set_ylim(top=ylim)
            plt.yticks(np.arange(0, (ylim+ylim/5), ylim/5), zorder=6)
            break
        else:
            ylim = maxy
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)     
    # EPE: Substituting the axis label with a txt label otherwise
    # the lable is below the grid (probably a bug). 
    # see my post on stackexchange: 
    # https://stackoverflow.com/questions/46242088/axis-label-hidden-by-axis-in-plot

    for it in np.arange(len(theta)):
        # first letter of the sector:
        ax.text(theta[it], ylim*1.15, index[it][0], va='center',
                ha='center', fontsize=ftx)
    ax.set_xticklabels('')
    ax.tick_params(axis='y', labelsize=ftx-2, zorder=4)
    ax.text(2,  ylim*1.3, titles_dct[prec], va='center',
                    ha='center', fontsize=(ftx+2), weight='bold')
    plt.close()
#    plt.show()
    return fig


def haversine_vec(lon1,lat1,lon_vec2,lat_vec2):
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
    # calculate haversine
    dlat = np.radians(lat_vec2) - np.radians(lat1)
    dlon = np.radians(lon_vec2) - np.radians(lon1)
    d = np.sin(dlat * 0.5) ** 2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat_vec2)) * np.sin(dlon * 0.5) ** 2
    h = 2 * AVG_EARTH_RADIUS * np.arcsin(np.sqrt(d))
    return h # in kilometers
#

def read_nuts_area(filenuts, calcall=False, nullnut=False, nutsall=None):
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
        WARNING by EPE: this function does not work andymore as originally 
        intended because the gridintersect structure has changed
        for example there is no nullnut anymore which has the same name for all 
        cells! the 'rect' par has not been tested.
    REFERENCES
    
    '''   
    nuts_info_all={}
    if(filenuts != 'rect'):
        nuts_def= filenuts +'.txt'
        nuts_info = pd.read_csv(nuts_def,delimiter="\t")
        nuts_info=nuts_info.dropna(axis=1,how='all')
        nutsnames=list(nuts_info.columns[~nuts_info.columns.isin(['POP','COL','ROW','AREA_km2','LAT','LON','CENTROID_X', 'CENTROID_Y', 'PERCENTAGE', 'POPULATION'])])
        #optional 'nut' comprising all grid points
        if calcall :
        #nutsnames.insert(0, 'ALL')
            nutsnames.insert(0, 'ALL_'+nutsnames[0])
            nuts_info[nutsnames[0]]=nutsnames[0]
        nuts_info['grid']=['_'.join(str(i) for i in z) for z in zip(nuts_info['COL'],nuts_info['ROW'])]
        if 'AREA_km2' in nuts_info.columns:
            nuts_area=pd.concat(map(lambda p: nuts_info['AREA_km2'],nutsnames),axis=1)
            #nuts_area.index=nuts_info['grid']
            nuts_area.columns=nutsnames
           #nuts_info=nuts_info[nutsnames]
        else:
            sys.exit("missing infos on grid cells area per nut")

        #aggregate data for each nut, create a dictionary
        nut_info_nut={}
        nut_info={}
        nullnut_dct={'NUTS_Lv0':'0', # at this level there is no background
                     'NUTS_Lv1':'1', 
                     'NUTS_Lv2':'11',
                     'NUTS_Lv3':'111'}
        for nut in nutsnames:
#            nut = 'NUTS_Lv0'
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
            if nullnut is True:
                for nutkey in set(nut_info_nut.index.get_level_values(0)):
                    if nutkey[2:]==nullnut_dct[nut]:
#                        print(nutkey)
                        nut_info_nut=nut_info_nut.drop(nutkey, level='nutname')
            nuts_info_all[nut]=nut_info_nut
       
    else:
        nuts_rect=nutsall
        nuts_rect.index=nuts_rect.index.droplevel(level=0)
        grid_inrect=nuts_rect.index
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']>=rect_coord['ll']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']<=rect_coord['ur']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']>=rect_coord['ll']['lat']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']<=rect_coord['ur']['lat']]
        nuts_rect=nuts_rect.loc[list(grid_inrect)]
        nuts_rect['nutname'] = 'rect'
        nuts_rect.set_index('nutname', append=True, inplace=True)
        nuts_info_all['rect']=nuts_rect.swaplevel(i=-2, j=-1, axis=0)
    return nuts_info_all

def find_target_info(areaid,areacells,x_y):
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


def dc_areasum(dc_all,narea,liv=1):
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
    dc_tot=pd.concat(list(map(lambda areaname:(dc_all[narea[areaname].index]*narea[areaname].values).sum(axis=1) ,narea.keys())),axis=1)
    dc_tot.columns=list(narea.keys())
    #overall sum including sum over snaps
    if isinstance(dc_all.index, pd.core.index.MultiIndex):
        alldc_snaps=dc_tot.groupby(level=[liv]).sum()
    else:
        alldc_snaps=dc_tot
    return alldc_snaps


def dc_snapaggregate(alldc_snaps,aggr_src=True):
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
    # aggregate sources as required in e-rep
    if(aggr_src==True):
        sources=pd.Series(['Transport','Industry','Industry','Industry','Agriculture','Residential','Other','Other','Other','Other'], #  Modified by Ema
                            index=['Nsnaps07','Nsnaps01','Nsnaps03','Nsnaps04','Nsnaps10','Nsnaps02','Nsnaps08','Nsnaps05','Nsnaps06','Nsnaps09'])
    else:
        sources=pd.Series(['Transport','Energy','Industry','Production','Agriculture','Residential','Offroad','Extraction','Industry','Waste'],
                            index=['Nsnaps07','Nsnaps01','Nsnaps03','Nsnaps04','Nsnaps10','Nsnaps02','Nsnaps08','Nsnaps05','Nsnaps06','Nsnaps09'])

    alldc=alldc_snaps
    alldc['source']=sources.loc[alldc_snaps.index]
    alldc.set_index('source', append=True, inplace=True)
    alldc=alldc_snaps.groupby(level=1).sum()
    if(aggr_src==True):  # EPE: changed order to have other at the end
                         # for better display
        alldc = alldc.reindex(['Transport', 'Industry', 'Agriculture', 'Residential',
                               'Other'])
    return alldc


def dc_increments(alldc,aggr_zones='fua'):   
    '''
    NAME
        calculate increments
    PURPOSE
        calculate increments
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
        Emanuela Peduzzi
    REFERENCES
    
    '''
    alldc_inc={}
    if(aggr_zones!='rect'):
        alldc_inc['ALL_NUTS_Lv0']=alldc['ALL_NUTS_Lv0']-alldc['NUTS_Lv0']
        if aggr_zones=='fua':
            if 'FUA_CODE' in alldc.columns:
                alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']-alldc['FUA_CODE']
                if 'CCITY_CODE' in alldc.columns:
                    # EPE: check if there is a commuting zone 
                    # some cities, like Liverpool don't have one
                    cond = all(alldc['FUA_CODE'].values==alldc['CCITY_CODE'].values)
                    if cond == True: 
                        # there is no commuting zone
                        alldc_inc['CCITY_CODE']=alldc['CCITY_CODE']
                    else: 
                        # there is a commuting zone
                        alldc_inc['FUA_CODE']=alldc['FUA_CODE']-alldc['CCITY_CODE']
                        alldc_inc['CCITY_CODE']=alldc['CCITY_CODE']
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
        # EPE: added to aggregate for FUAs only (and used also when the area is too small)
        if aggr_zones=='fuaonly':
            if 'FUA_CODE' in alldc.columns:
                alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']-alldc['FUA_CODE']
                alldc_inc['FUA_CODE']=alldc['FUA_CODE']
            else:
               alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']

    else:
        alldc_inc['ALL_NUTS_Lv0']=alldc['ALL_NUTS_Lv0']-alldc['rect']
        alldc_inc['rect']=alldc['rect']

    alldc_inc=pd.DataFrame.from_dict(alldc_inc)
    return alldc_inc

def sherpa_model(pr,model_idx,alldists,emissiondelta,inner_radius=False):
    '''
    NAME
        Core of sherpa model calculation
    PURPOSE
        Calculate the delta conc for one precursor on all snaps uding a predefined set of grid point for emission reduction (area)
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
    
    REFERENCES
    
    '''
    emidelta=emissiondelta.loc[pr]
    #fastest
    window_all=(1.+alldists)**(-model_idx.loc['omega'][pr])
    dc_cells=model_idx.loc['alpha'][pr]*(emidelta*window_all)
    if inner_radius is not False : #if not zero then use info on flatweight
        outer_cells=emidelta.sum()[(alldists>inner_radius) & (emidelta.sum()>0)].index 
        if len(outer_cells)>0:
       # print ("fo prec",pr,"applying flatweight on",len(outer_cells),"points in area, with",len(alldists[alldists>inner_radius].index),
       # "cells beyond inner radius but",len(alldists[emidelta.sum()==0].index),"zero delta emissions")
            weighted_emissions_flat = model_idx.loc['alpha'][pr]*(model_idx.loc['flatWeight'][pr] * emidelta[outer_cells].sum(axis=1))
    #disaggregate the value on all involved points
            flat_value=weighted_emissions_flat/len(outer_cells)
            flat_matrix=pd.DataFrame(index=outer_cells, columns=flat_value.index)
            flat_matrix=flat_matrix.fillna(flat_value)
            dc_cells.loc[:,outer_cells]=flat_matrix.transpose()
    return dc_cells



def df2mat(dfdata):
    '''
    NAME
        Convert a pd.seroes in a 2D matrix
    PURPOSE
        Convert a pd.seroes in a 2D matrix
    PROGRAMMER(S)
        Denise Pernigotti starting from Bart routine for writing delta conc netcdf in module 1
    REVISION HISTORY
    
    REFERENCES

    '''
    grids=pd.DataFrame(dfdata.index.str.split('_',expand=True).tolist(), columns=['x','y'], index=dfdata.index)
    grids['x']=pd.to_numeric(grids['x'])
    grids['y']=pd.to_numeric(grids['y'])
    grids['var']=dfdata
    dfmat=grids.pivot(index='y', columns='x', values='var')
    dfmat=dfmat.as_matrix()
    return dfmat

def module7(emissions_nc, concentration_nc, natural_dir, model_nc, fua_intersect_dir, nuts_intersect_dir, dbf_dir, targets_txt, outdir, aggr_zones_in, pollutant, outfig='png', normalize=True):
    '''
    NAME
        Module 7
    PURPOSE
        For each target in a list loaded from file target.txt produce data and plots for e-reporting:
            --> calculate the variation of the concentration dc due to each grid point 100% emission reduction of each precursor;
            --> the dc increments on user defined areas (nuts/fua etc.) in calculated;
            --> a png plot with the concentration increments flattened on two the smaller areas is produced
            --> these data are then aggregated  for snaps and precursors producing 2 html tables;
    PROGRAMMER(S)
        Denise Pernigotti June 2017, Peduzzi Emanuela (EPE) January 2018,
    REVISION HISTORY
        From last version of module7_custom.py
    
    REFERENCES
    
    '''
    # EPE rect_coord not used for the moment
#    rect_coord=None
      
    ################ default user definition
    aggr_src=True #set to True for aggregatin sources as in the atlas
    # EPE: if pollutant is NOx there are no corresponding natural concentrations files
    if pollutant == 'NOx':
        include_natural=False
        poll_name = 'NO2'
    else: 
        include_natural=True #include or not natiral PM concentration
        poll_name = pollutant
    aggr_sd=True # aggregate salt and dust in natural (EPE)
    print_areainfo=False #if true check and print which grid points in fua/gcity are actually not modelled
    #############
    
    # Cehck consistency of pollutant and input files---------------------------
    if poll_name in emissions_nc and concentration_nc and natural_dir and model_nc: 
        print('OK: Pollutant name and input files are consistent')
    else: 
        sys.exit('ERROR: pollutant name and input files are not consistend')
        
#    rect_txt = intersect_dir+'/gridint_rect' EPE: not sure what this is for now
    # EPE new grid intersect:
    grd_fua_txt = fua_intersect_dir + 'grid_intersect' 
    grd_nuts_txt = nuts_intersect_dir +'grid_intersect'

    #list of true names for areas IDs    
    codes_names=pd.Series({'NUTS_Lv0':'fua/FUA_2013_WGS84_Lv0.dbf', # 
                           'NUTS_Lv1':'nuts/NUTS_2013_WGS84_Lv1.dbf', # 
                           'NUTS_Lv2':'nuts/NUTS_2013_WGS84_Lv2.dbf',
                           'NUTS_Lv3':'nuts/NUTS_2013_WGS84_Lv3.dbf',
                           'FUA_CODE': 'fua/FUA_2013_WGS84_Lv2.dbf',
                           'CCITY_CODE':'fua/FUA_2013_WGS84_Lv3.dbf'}) # as in Core City code

    codes_txt={k: dbf_dir + codes_names[k] for k in codes_names.keys()}

    #list of points used as receptors. Only 'id','lon','lat' and 'Station name' are required columns
    ###############################################

    #read netcdf files, put the data in multindex dataframes and check consistency
    emissions = read_nc(emissions_nc)

    concentration = read_nc(concentration_nc)
    model= read_nc(model_nc)

    if include_natural:
        pmsize=pollutant[-2:]
        salt_nc = natural_dir + 'pSALT-'+pmsize+'-basecase.nc'  
        dust_nc = natural_dir + 'pDUST-'+pmsize+'-basecase.nc'  
        salt= read_nc(salt_nc)
        dust= read_nc(dust_nc)

    #check consistency of model and emissions
    if model.loc['coord'].astype(np.float32).equals(emissions.loc['coord'].astype(np.float32)):
        print ('OK: latitude and longitude in matrices emissions and model are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")

    #check consistency of model and concentrations
    if model.loc['coord'].astype(np.float32).equals(concentration.loc['coord'].astype(np.float32)):
    #if model.loc['coord'].equals(concentration.loc['coord']):
        print ('OK: latitude and longitude in matrices model and conc are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")

    #get inner radius and check if the flat weight option is activated
    for nameatt in Dataset(model_nc,'r').ncattrs():
        if nameatt[0:6] == 'Radius':
            radiusofinfluence = nameatt
    inner_radius = int(getattr(Dataset(model_nc,'r'), radiusofinfluence))
    if 'flatWeight' in model.index.levels[0]:
        print('Flatweight approximation is ON with Radius of influence of '+str(inner_radius) +' grid cells')
#        if inner_radius >=200: # @todo EPE removed this condition, this needs to be checked (EPI?)
#            sys.exit("Something wrong, radius of influence is >=200 grid cells")
    else:
#        if inner_radius is not 200:
#            sys.exit("Something wrong, flatweight approximation is deactivated but inner_radius is not 200 grid cells")
        print('Flatweight approximation is OFF')
        inner_radius=False

    #save coord info
    coordinates=emissions.loc['coord',].transpose()
    lat_vec=np.unique(coordinates['lat'])
    lon_vec=np.unique(coordinates['lon'])
    #find grid resolution in degrees
#    dlon_res=min(lon_vec[1:]-lon_vec[:-1])
#    dlat_res=min(lat_vec[1:]-lat_vec[:-1])

    # grab nuts/areas info from txt
    nuts_info=read_nuts_area(grd_nuts_txt,calcall=True)    
    nuts_info_fuas=read_nuts_area(grd_fua_txt, nullnut=True) 
    nuts_info_fuas['CCITY_CODE']=nuts_info_fuas['NUTS_Lv3']
    nuts_info_fuas['FUA_CODE']=nuts_info_fuas['NUTS_Lv2']
    del nuts_info_fuas['NUTS_Lv3']
    del nuts_info_fuas['NUTS_Lv2']
    del nuts_info_fuas['NUTS_Lv1']
    del nuts_info_fuas['NUTS_Lv0']  
    nuts_info.update(nuts_info_fuas)
    # EPE, only one instance of read_nuts_area should have calcall=True, original 
    # version allowed to set nullnut='LAND000', but now this nut name does not exist 
    # anymore.
    #, 
#    if(aggr_zones_in=='rect'):
#        nuts_info.update(read_nuts_area('rect',nutsall=nuts_info['ALL_NUTS_Lv0'].copy()))
##    nuts_info.keys()
#
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
    if receptors['lon'].min()<lon_vec.min() or receptors['lon'].max()>lon_vec.max():
        sys.exit("some receptors have longitude outside the domain")
    if receptors['lat'].min()<lat_vec.min() or receptors['lat'].max()>lat_vec.max():
        sys.exit("some receptors have latitude outside the domain")

    #grab true names for each area (Need the codec for the NUTS)
    area_names_long={k: (Dbf5(codes_txt[k], codec='latin1').to_dataframe()[['NUTS_ID','NAME_ASCI']].set_index('NUTS_ID')) for k in codes_names.keys()} 
  
    #reduce string length if needed
    area_names={}
    countries=np.unique(nuts_info['NUTS_Lv0'].index.get_level_values(0))
    for k in area_names_long.keys():
        area_names[k]=area_names_long[k]['NAME_ASCI'].apply(name_short)
    area_names['ALL_NUTS_Lv0']=pd.Series({'ALL_NUTS_Lv0' : 'Europe', 'other' : 'other'},name='EU')
    if(aggr_zones_in=='rect'):
        area_names['rect']=pd.Series({'rect' : 'Domain', 'other' : 'other'},name='EU')
    totalname='Total'
#    ########
    #optimize calculations removing unused data
    model=model.drop('coord',level=0)
    emissions=emissions.drop('coord',level=0)
    #fake sum to remove one level
    concentration=concentration.groupby(level=[0]).sum()
    concentration=concentration.loc['conc']

    #remove not modelled points (SR was build only using reductons in points in grid intersect)
    modelled_emissions_idx=list(set(emissions.columns).intersection(nuts_info['NUTS_Lv0'].index.get_level_values(1)))
    emissions=emissions[modelled_emissions_idx]
    print('in emissions keeping ',len(modelled_emissions_idx), ' grid points used in SR model training')
    print('the considered emissions are from '+str(len(countries))+' european or near european countries')
    print(countries)
    #check coherence with nuts_info grid points (must be all and only the ones in emissions)
    info_grids={}
    for area in nuts_info.keys():
        index_diff=set(nuts_info[area].index.get_level_values(1)).difference(modelled_emissions_idx)
        if len(index_diff)>0:
            idx = pd.IndexSlice
            print('There are '+str(len(index_diff)) +
            ' points in area '+area+' which are not modelled an will be removed from area summing up')
            if print_areainfo:
                info_grids_area=coordinates.loc[index_diff]
                info_grids_area=info_grids_area.merge(nuts_info[area].reset_index(level=0),how='left',left_index=True, right_index=True)
                info_grids[area]=info_grids_area.merge(area_names[area].reset_index(level=0),how='left',left_on='nutname',right_on=0,left_index=True)
                info_grids[area]=info_grids[area].set_index(info_grids_area.index)
                info_grids[area]= info_grids[area].drop([0], 1)
                print(info_grids[area])
            nuts_info[area]=nuts_info[area].loc[idx[:,modelled_emissions_idx],:]


    #remove not modelled grid points
    model_sum=model.loc['alpha'].sum()
    model_points=model_sum[model_sum>0].index
    print ("in model removing ",len(model_sum[model_sum==0].index)," grid points where alpha is zero")
    model=model[model_points]

    #remove not modelled precursors
    model_sumgrid=model.loc['alpha'].sum(axis=1)
    precursors=list(model_sumgrid[model_sumgrid>0].index)
    print ("in model keep ",list(precursors)," modelled precursors")
    model=model.loc[(slice(None),precursors),:]
    emissions=emissions.loc[(precursors,slice(None)),:]
    #get measurement unit for emissions
    emi_units=list(map(lambda p: getattr(Dataset(emissions_nc,'r').variables[p], 'units'),precursors))
    if len(set(emi_units))!=1:
        sys.exit("not all precursors have the same unit mass in emissions")

    alldist_km=pd.concat(list(map(lambda st: haversine_vec(receptors.loc[st,'lon'],receptors.loc[st,'lat'],coordinates['lon'],coordinates['lat']),receptors.index)),axis=1)
    receptors['id']=receptors['id'].str.strip()
    receptors['target_idx']=alldist_km.idxmin()
    receptors['dist_km']=alldist_km.min()
    receptors['lon_grid']=pd.Series({st:coordinates.loc[receptors.loc[st,'target_idx'],'lon'] for st in receptors.index})
    receptors['lat_grid']=pd.Series({st:coordinates.loc[receptors.loc[st,'target_idx'],'lat'] for st in receptors.index})
    receptors['model_conc']=concentration[receptors['target_idx']].values
    print("the selected targets correspond to modelled grid points x_y ")
    print(receptors['target_idx'].values)
    count_idx=receptors.pivot(columns='target_idx', values='target_idx').count()
    if count_idx.max()>1:
        print('There are duplicates in target_idx:' + ', '.join(list(count_idx.loc[count_idx>1,].index)))
    a =count_idx.reset_index()
    a.rename(columns={0: 'duplicates'}, inplace=True)
    receptors=pd.merge(receptors,a)
    receptors.index=receptors['id']
    receptors.drop('id', axis=1, inplace=True)

    dc_inc_all={}
    dc={}
    target_allinfo={}
    #calculate diftsnces first in order to save calculation time, targets nedd to be in  memory limit of the machine (say at most about 10000)
    dists_array=distance.cdist(coordinates.loc[count_idx.index,['x','y']],coordinates.loc[emissions.columns,['x','y']], metric='euclidean')

    for ix,idx in enumerate(count_idx.index):
        #For the selected point find all information on area  ids
        target_info=pd.concat(list(map(lambda areaid: find_target_info(areaid,nuts_info,idx),nuts_info.keys())),axis=1).transpose()
        target_info['areaname']=pd.Series(dict(zip(target_info.index,map(lambda x: ''.join(list(area_names[x].loc[target_info.loc[x,'areaid']])).title(),target_info.index))))
        #select grid points in the same area as the target, with their percentual area
        narea=pd.concat(list(map(lambda areaname:nuts_info[areaname].loc[target_info.loc[areaname,'areaid']]['parea'],target_info.index)),axis=1)
        narea.columns=target_info.index
        narea=narea.fillna(0)
        #narea = dict(zip(target_info.index, narea))

        #Core of the sherpa calculations: calculate dc due to each gridpoint for each precursor and macrosector
        dc[idx]=pd.concat(list(map(lambda p: sherpa_model(p, model[idx],dists_array[ix,:],emissions,inner_radius),precursors)))
        dc[idx].index=emissions.index

        #aggregate emissions per precursor and area
        nareasurf=pd.concat(list(map(lambda areaname:nuts_info[areaname]['area'].loc[target_info.loc[areaname,'areaid']],target_info.index)),axis=1)
        nareasurf.columns=target_info.index
        nareasurf=nareasurf.fillna(0)
        
        # EPE: change aggregations in case of exceptions:  
        # check the area of the citycore, 
        # if it is less then 300 the aggregation zone is switched to 
        # fua_only 
        print('checking aggregation zone', aggr_zones_in, 'for', idx)
        if aggr_zones_in == 'fua' or aggr_zones_in == 'fuaonly':
            aggr_zones = aggr_zones_in
            if 'CCITY_CODE' not in list(target_info.index) or 'FUA_CODE' not in list(target_info.index):
                aggr_zones='nuts'
                print('The selected point is not in a FUA, changing aggregation to NUTS')
            if 'CCITY_CODE' in list(target_info.index):
                print('the area of the city core is:', nareasurf['CCITY_CODE'].sum()) 
                if nareasurf['CCITY_CODE'].sum() <= 300:
                    print('The city core is', nareasurf['CCITY_CODE'].sum(), '<= 300 km2, therefore the whole FUA is considered')
                    if aggr_zones_in == 'fua':
                        aggr_zones ='fuaonly'         
        else: 
            aggr_zones = aggr_zones_in
            
        print('using aggragation zones as ', aggr_zones)
            #define the aggregation and increments calculation type depending on aggr_zones
        if aggr_zones=='fua':
            wantedorder=pd.Series( {'3' : 'ALL_NUTS_Lv0', '2' : 'NUTS_Lv0', '1' :'FUA_CODE','0':'CCITY_CODE'})
        elif aggr_zones=='fuaonly':
            wantedorder=pd.Series( {'1' : 'FUA_CODE', '2':'NUTS_Lv0', '3':'ALL_NUTS_Lv0'})
        elif aggr_zones=='nuts':
            wantedorder=pd.Series( {'0' : 'NUTS_Lv3','1' : 'NUTS_Lv2', '2' : 'NUTS_Lv1', '3' :'NUTS_Lv0','4':'ALL_NUTS_Lv0'})
        elif aggr_zones=='rect':
            wantedorder=pd.Series( { '1' :'rect','0':'ALL_NUTS_Lv0'})

        
        emi_sum=pd.concat(list(map(lambda p: dc_snapaggregate(dc_areasum(emissions.loc[p],nareasurf),aggr_src),precursors)),axis=0,keys=precursors)
        emi_inc=dc_increments(emi_sum,aggr_zones)
        emi_inc[totalname]= emi_inc.sum(axis=1)

        #aggregate dc per precursor, area, calculate increments and relative values
        alldc=dc_snapaggregate(dc_areasum(dc[idx],narea),aggr_src)
        dc_inc=dc_increments(alldc,aggr_zones)*100./concentration[idx]
        area_present=dc_inc.columns
        dc_inc[totalname]= dc_inc.sum(axis=1)

        #aggregate results per precursor
        alldc_prec=dc_areasum(dc[idx],narea,liv=0)
        dc_inc_p=dc_increments(alldc_prec,aggr_zones)*100./concentration[idx]
        dc_inc_p[totalname]= dc_inc_p.sum(axis=1)

        wantedorder_present=pd.Series(list(filter(lambda x: x in area_present, wantedorder))).to_frame(name='areaid')
        wantedorder_present['areaname']=pd.Series(list(target_info.loc[wantedorder_present['areaid'],'areaname']))
         #check for duplicated names in nuts
        if 'NUTS_Lv1' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv1','areaname']='1_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv1','areaname']
        if 'NUTS_Lv2' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv2','areaname']='2_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='NUTS_Lv2','areaname']

         #rename fua elements (same names as city core)
        if 'FUA_CODE' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']= wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']
            # wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']='fua_'+wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']

        #build a dataframe with constant values in the smallest areas (excluding 'ALL_NUTS_Lv0' and 'NUTS_Lv0')
        smallareas=list(set(area_present)-set(['ALL_NUTS_Lv0','NUTS_Lv0']))

        #avoid double counting of grid increments
        narea_inc=dc_increments(narea,aggr_zones)

        if len(smallareas)>0:
        #make a copy with 0 and 1 only for plotting
            narea_bin=narea[smallareas].copy()
            narea_bin=narea_bin.mask(narea_bin>0.5,1)
            narea_bin=narea_bin.mask(narea_bin<=0.5,0)
            narea_id=narea_bin.sum(axis=1)

            #dc_inc_flat= pd.concat(list(map(lambda p: narea_inc[p]*(dc_inc[p].sum()),smallareas)),axis=1).sum(axis=1)
            #dc_inc_flat= dc_inc_flat.reindex(index=dc[idx].columns)

        #set up appropriate names for the areas
        wantedorder_present=wantedorder_present.append(pd.Series({'areaid':totalname, 'areaname':totalname}), ignore_index=True)
        dc_inc=dc_inc[wantedorder_present['areaid']]
        #add natural sources
        if include_natural:
            # EPE: modified to allow aggregation of natural contribution
            if aggr_sd == False: 
                natural=pd.DataFrame(0, index=['Salt','Dust'], columns=dc_inc.columns)
                natural.loc['Dust',totalname]=dust.loc['pDUST-'+pmsize,idx].values*100./concentration[idx]
                natural.loc['Salt',totalname]=salt.loc['pSALT-'+pmsize,idx].values*100./concentration[idx]
            elif aggr_sd == True: 
                natural=pd.DataFrame(0, index=['Natural'], columns=dc_inc.columns)
                natural.loc['Natural', totalname] = (
                        dust.loc['pDUST-'+pmsize,idx].values
                        + salt.loc['pSALT-'+pmsize,idx].values) * 100./concentration[idx]           
            dc_inc=dc_inc.append(natural)

        #plots
        fig={}
        if aggr_zones=='fua' or  aggr_zones=='fuaonly':
            az_name = 'fua'
        else: 
            az_name = 'nuts'
        
        fig[1]=plot_bar(dc_inc,wantedorder_present,totalname, normalize=normalize)
        if aggr_zones=='fua' or  aggr_zones=='fuaonly':
            for ip,p in enumerate(precursors):
                fig[2+ip] = plot_polar(emi_sum, p, wantedorder_present)
        for ids in list(receptors[receptors['target_idx']==idx].index):
            fig[1].savefig(outdir+'\\'+ids+'_'+pollutant+'_'+az_name+'_conc.'+ outfig, dpi=1000, bbox_inches='tight')

            if len(smallareas)>0:
                if aggr_zones=='fua' or  aggr_zones=='fuaonly':
                    for ip,p in enumerate(precursors):
                        if fig[2+ip]:  
                            fig[2+ip].savefig((outdir+'\\'+ids+'_'+p+'_emi.'+outfig), dpi=300, bbox_inches='tight')
            dc_inc.to_csv(outdir+'\\'+ids+'_'+pollutant+'_'+az_name+'_total_table.csv')
            dc_inc_all[ids]=dc_inc.transpose()
            target_allinfo[ids]=target_info.transpose()
            # EPE: make legend and save it
            figleg = plt.figure(figsize=figsize(0.1))#
            axleg = figleg.add_subplot()
            axleg = plt.subplot()
            axleg.set_axis_off()
            ax1 = fig[1].get_axes()
            handles1, labels1 = ax1[0].get_legend_handles_labels()
            # add first letter to labels of sectors
            for labit in np.arange(len(labels1)):
                labels1[labit]= labels1[labit][0] +' - '+labels1[labit]
                
            if aggr_zones=='fua' or  aggr_zones=='fuaonly':              
                ax2 = fig[2].get_axes()
                handles2 = ax2[0].lines
                labels2 = [handles.get_label() for handles in handles2]                       
                dct_names = {'CCITY_CODE': 'City',
                             'FUA_CODE': 'Greater city'} 
                newlabels = [dct_names[label] for label in labels2]
            else: 
                handles2=[]
                newlabels=[]
                
            axleg.legend(handles=handles1+handles2, labels=labels1+newlabels, fontsize=10, loc = 'center', frameon=False) 
            figleg.savefig((outdir+'\\'+ids+'_'+pollutant+'_'+az_name+'_'+'_legend.'+outfig), dpi=300, bbox_inches='tight')
            plt.close('all')

    #summarize info on grid points
    reform = {(outerKey, innerKey): values for outerKey, innerDict in target_allinfo.items() for innerKey, values in innerDict.items()}
    target_allinfo=pd.DataFrame(reform).transpose()
    a =target_allinfo.reset_index()
    a.rename(columns={'level_0': 'id', 'level_1': 'area'}, inplace=True)
    b =receptors[['station name','target_idx','duplicates','model_conc','lon','lon_grid','lat','lat_grid','dist_km']].reset_index()
    summary= pd.merge(b, a)
    summary = summary.set_index(['id','area'])
    summary.to_csv(outdir+'\\AAA_summary_info.csv')

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

##-----------------------------------------------------------------------------    
#    # Extract data for the Covenant of Mayors analysis EPE
#    ## produce one file for each allocation level
#    df_res = {}
#    cities_com_coord = pd.read_csv(('D:/sherpa.git/Sherpa/com20170706/' + 'Cities_CoM_coord.csv'), index_col=[2])
#    cities_com_coord.index.rename('id', inplace = True)
#    del cities_com_coord['oid']
#    del cities_com_coord['seap_id']
#    if aggr_zones=='fua' or aggr_zones=='fuaonly':
#        levels = ['CCITY_CODE', 'FUA_CODE', 'NUTS_Lv3', 'NUTS_Lv2',
#                       'NUTS_Lv1', 'NUTS_Lv0', 'ALL_NUTS_Lv0', 'Total']
#    else: 
#        levels = ['NUTS_Lv3', 'NUTS_Lv2',
#                      'NUTS_Lv1', 'NUTS_Lv0', 'ALL_NUTS_Lv0', 'Total']
#        
#    for level in levels:        
#            df_res[level]  = pd.DataFrame(columns=set(dc_inc_all.index.get_level_values(1)))
#
#    df_res_tot=pd.DataFrame(0, index=set(dc_inc_all.index.get_level_values(0)), columns=set(dc_inc_all.index.get_level_values(1)))
#    for indx in set(dc_inc_all.index.get_level_values(0)):
#        for level in levels:    
#            if not np.isnan(dc_inc_all.loc[indx][level]).all():
#                df_res_tot.loc[indx]=df_res_tot.loc[indx]+dc_inc_all.loc[indx][level] 
#                df_res[level].loc[indx]=df_res_tot.loc[indx]
#    for level in levels:        
#        results = cities_com_coord.join(receptors.join(df_res[level]))              
#        results.to_csv(outdir + 'sourceall{}_{}.csv'.format(pollutant, level))
##-----------------------------------------------------------------------------

if __name__ == '__main__':

    pass
