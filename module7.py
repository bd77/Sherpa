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
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np
import pandas as pd
def read_nc(nc_file):
    #nc_file='input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    #nc_file='input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    nc_data = Dataset(nc_file, 'r')
    nc_dims = [dim for dim in nc_data.dimensions]
    nc_vars = [var for var in nc_data.variables]
    #sometimes the latitude is written just with lat as in model data
    latname=list(filter(lambda x: x in nc_vars, ['Lat','lat','latitude']))[0]
    #latname=list(filter (lambda x: 'lat' in x, nc_vars))[0]
    lats = nc_data.variables[latname][:]
    lonname=list(filter(lambda x: x in nc_vars, ['Lon','lon','longitude']))[0]
    #lonname=list(filter (lambda x: 'lon' in x, nc_vars))[0]
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
            #nznames=[ncz + s for s in map(str,range(1,len(nc_data.dimensions[ncz])+1))]
    #create an index with lat and lon
    #latrep=map(str, np.repeat(lats,len(lons)))
    #lonrep=map(str, np.tile(lons,len(lats)))
    #trasform variables arrays in vectors
    #allvar={'lat_lon':map(lambda (x,y): x+'_'+y, zip(latrep, lonrep))}
    #create lat and lon info
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
        #allvarnc[var]=allvarnc[var].transpose()
        #index_var = pd.MultiIndex.from_tuples(zip(np.repeat(var,len(nz)),nznames), names=['vars', ncz])
        #allvar[var].columns=index_var
    reform = {(outerKey, innerKey): values for outerKey, innerDict in allvar.items() for innerKey, values in innerDict.items()}
    df=pd.DataFrame(reform)
    return df.transpose()

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
    Write a 2D matrix or a dictionary of 2D matrices on a netcdf file
PURPOSE
    CWrite a 2D matrix or a dictionary of 2D matrices on a netcdf file
PROGRAMMER(S)
    Denise Pernigotti starting from Bart routine for writing delta conc netcdf in module 1
REVISION HISTORY

REFERENCES

'''
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
def write_dict_nc(dc_dic,lats,lons,unit,filenc):
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

'''
NAME
    Spatial plot a column in a panda df
PURPOSE
    Spatial plot a column in a panda df, winfos covers grid_grads degrees
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY

REFERENCES

'''
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#from pylab import *
def plot_dict(dc_dic,idx,coord_emissions,grid_grads=3,unit='%'):
    lat_emissions=np.unique(coord_emissions['lat'])
    lon_emissions=np.unique(coord_emissions['lon'])
    #grid_grads define the area for the plot
    #x=int(idx.split('_')[0])-1
    #y=int(idx.split('_')[1])-1
    lat_0=coord_emissions.loc[idx,'lat']
    lon_0=coord_emissions.loc[idx,'lon']
    ll_lon=lon_0-grid_grads
    ll_lat=lat_0-grid_grads
    ur_lon=lon_0+grid_grads
    ur_lat=lat_0+grid_grads
    lon, lat = np.meshgrid(lon_emissions, lat_emissions)
    dfmat=df2mat(dc_dic)
    plt.close('all')
    fig=plt.figure()
    m = Basemap(resolution='i',projection='cyl',
            llcrnrlon=ll_lon, llcrnrlat=ll_lat,
            urcrnrlon=ur_lon, urcrnrlat=ur_lat,
    #        llcrnrlon=min(lons), llcrnrlat=min(lats),
    #        urcrnrlon=max(lons), urcrnrlat=max(lats),
            area_thresh = 0.1,lat_0=lat_0,lon_0=lon_0)

    m.drawcoastlines();
    #m.drawstates()
    m.drawcountries();
    # draw filled contours.
    x, y = m(lon, lat)
    #clevs = [5,10,15,20,25,30,35,40,45,50,55,60]
    #clevs=[0.9,1.1,1.9,2.1,2.9,3.1,3.9,4.1]
    clevs=[0,1,2,3,4]
    #cs = m.contour(x,y,dfmat,clevs);
    cs = m.contour(x,y,dfmat,clevs,linewidths=1,colors=('r', 'green', 'blue','grey','black'));
    # add colorbar.
    #cbar = m.colorbar(cs,location='bottom',pad="5%")
    #cbar.set_label(unit)
    #drow target
    x,y = m(lon_0, lat_0)
    m.plot(x, y, 'x', markersize=6,color='k',fillstyle='none')
    #plt.title(idx)
    #plt.savefig(fileout)
    #plt.close(fig)
    #plt.show()
    return fig


from matplotlib.ticker import AutoMinorLocator # EMA this should be at the top


def plot_bar(dfdata, varplot, totalname, plot_opt='perc',
             x_label='% of total mass', leg_loc='lower right', ftx=16):
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
        Revised by @author: peduzem 20/09/2017
    REFERENCES
    '''

    dfdata = dfdata[varplot['areaid']].transpose()
    varnames = dfdata.columns
    areanames = varplot['areaname']
    areanames.index = varplot['areaid']
    varplot_left = pd.Series(np.zeros(len(dfdata.index)), index=dfdata.index)
    addsum = dfdata.sum(axis=1)

    for ind in range(1, len(dfdata.index)):
        varplot_left[ind] = varplot_left[ind-1]+addsum[dfdata.index[ind-1]]
    if totalname in varplot['areaid'].values:
        varplot_left[totalname] = 0.

    # colors by EPI
    colors = {'PPM': 'blue', 'SOx': 'gold', 'NOx': 'red', 'NH3': 'green',
              'NMVOC': '#black', 'Traffic': 'red', 'Energy': 'blue',
              'Industry': 'gold', 'Production': '#8B4789', 'Waste': '#00FFFF',
              'Agriculture': 'green', 'Residential': 'blue',
              'Offroad': '#8470ff', 'Extraction': '#00FF00',
              'Other': 'skyblue', 'Salt': '#ccffe5', 'Dust': '#ffffcc'}

#    colors={'PPM':'grey','SOx':'#ff0000','NOx': '#0000ff','NH3':'#b22222','NMVOC':'#8470ff',
#    'Traffic': '#0000ff','Energy':'#999999','Industry':'#ff0000','Production':'#8B4789','Waste':'#00FFFF','Agriculture':'#b22222','Residential':'#32cd32',
#    'Offroad':'#8470ff','Extraction':'#00FF00','Other':'#bebebe','Salt':'#ccffe5','Dust':'#ffffcc'}

    # Create the general blog and the "subplots" i.e. the bars
    plt.close('all')
    f, ax1 = plt.subplots(1, figsize=(12, 6))  # modified by EPI
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
             dfdata[varnames[0]],
             # set the width
             height=bar_width,
             # with pre_score and mid_score on the bottom
             left=list(varplot_left),
             label=varnames[0],
             # with alpha 0.5
             # alpha=0.5,
             # with color
             color=colors[varnames[0]])

    # Create a bar plot, in position bar_1
    for ivar in range(1, len(varnames)):
        ax1.barh(bar_l,
                 # using the mid_score data
                 dfdata[varnames[ivar]],
                 # set the width
                 height=bar_width,
                 # with pre_score and mid_score on the bottom
                 left=list(varplot_left+dfdata[
                           varnames[range(0, ivar)]].sum(axis=1)),
                 # with the label post score
                 label=varnames[ivar],
                 # with alpha 0.5
                 # alpha=0.5,
                 # with color
                 color=colors[varnames[ivar]])

    # set the y ticks with names
    plt.yticks(tick_pos, areanames, size=11, rotation=0)

    # set minor x ticks - added EMA
    minorLocator = AutoMinorLocator()
    ax1.xaxis.set_minor_locator(minorLocator)

    # Set the label and legends
    ax1.set_xlabel(x_label, fontsize=ftx)
    # ax1.set_ylabel("Areas")
    plt.legend(loc=leg_loc)

    # Set a buffer around the edge
    plt.ylim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
    return f


def figsize(scale):  # Added by Ema for spyder plot
    fig_width_pt = 390                              #
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size


def plot_spyder(emi_sum, prec):  # Added by Ema for spyder plot
    """Spyder plot for emissions, only if the point is in GCITY
    @author: peduzem

    """
#    prec = 'PPM'
    if 'GCITY_CODE' in emi_sum.columns:
        sect_display = emi_sum['GCITY_CODE'].loc[prec].index.get_level_values(0)

        index = [sect_display[1], sect_display[3], sect_display[4],
                 sect_display[0], sect_display[2]]
        df_city = emi_sum['GCITY_CODE'].loc[prec].reindex(index)/1000
        fig = plt.figure()
#        fig = plt.figure(figsize=figsize(1))
        ax = fig.add_subplot(111, projection="polar")
        ax.grid(True)
        ax.yaxis.grid(color='#dddddd')
        ax.xaxis.grid(color='#dddddd')
        colors = ['gold', 'tomato']
        for spine in ax.spines.values():
            spine.set_edgecolor('None')
        theta = np.arange(len(df_city))/float(len(df_city))*2.*np.pi
        ax.set_xticks(theta)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_horizontalalignment('center')
            tick.label1.set_verticalalignment('top')

        df = []

        if 'FUA_CODE' in emi_sum.columns:
            df_fua = emi_sum['FUA_CODE'].loc[prec].reindex(index)/1000
            df.append([df_fua])
        df.append([df_city])
        plots = []
        for it, dfdata in enumerate(df):
            line, = ax.plot(theta, df[it][0], color=colors[it],
                            marker="o", label=None, zorder=3)  # , zorder = -3)
            ax.fill(theta, df[it][0], color=colors[it], alpha=0.3,
                    zorder=3)
            plots.append(line,)

        def _closeline(line):
            x, y = line.get_data()
            x = np.concatenate((x, [x[0]]))
            y = np.concatenate((y, [y[0]]))
            line.set_data(x, y)
        [_closeline(l) for l in [plots[0], plots[1]]]

        ax.set_rlabel_position(90)  # get radial labels away from plotted line

        maxy = max(max(df_fua), max(df_city))
        for it in np.arange(len(theta)):
            ax.text(theta[it], maxy*1.15, index[it], va='center',
                    ha='center', fontsize=12)
        ax.set_xticklabels('')

        plt.legend(handles=[plots[0], plots[1]], labels=['FUA', 'city core'],
                   bbox_to_anchor=(0, 1), loc='center left')
        ax.tick_params(axis='y', labelsize=10, zorder=4)
        plt.title("{} emissions [kton/y]".format(prec), fontsize=14, y=1.2)
        return fig



##


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
import numpy as np
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
import pandas as pd
def read_nuts_area(filenuts,calcall=False,nullnut=None,nutsall=None):
    nuts_info_all={}
    if(filenuts != 'rect'):
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
def dc_snapaggregate(alldc_snaps,aggr_src=True):
    #aggregate sources as required in e-rep
    if(aggr_src==True):
        sources=pd.Series(['Traffic','Industry','Industry','Industry','Agriculture','Residential','Other','Other','Other','Other'], #  Modified by Ema
                            index=['Nsnaps07','Nsnaps01','Nsnaps03','Nsnaps04','Nsnaps10','Nsnaps02','Nsnaps08','Nsnaps05','Nsnaps06','Nsnaps09'])
    else:
        sources=pd.Series(['Traffic','Energy','Industry','Production','Agriculture','Residential','Offroad','Extraction','Industry','Waste'],
                            index=['Nsnaps07','Nsnaps01','Nsnaps03','Nsnaps04','Nsnaps10','Nsnaps02','Nsnaps08','Nsnaps05','Nsnaps06','Nsnaps09'])

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
    if(aggr_zones!='rect'):
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
### COM to remove start **
        if aggr_zones=='fuaonly':
            if 'FUA_CODE' in alldc.columns:
                alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']-alldc['FUA_CODE']
                alldc_inc['FUA_CODE']=alldc['FUA_CODE']
            else:
               alldc_inc['NUTS_Lv0']=alldc['NUTS_Lv0']
### COM to remove end **
    else:
        alldc_inc['ALL_NUTS_Lv0']=alldc['ALL_NUTS_Lv0']-alldc['rect']
        alldc_inc['rect']=alldc['rect']

    alldc_inc=pd.DataFrame.from_dict(alldc_inc)
    return alldc_inc

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
import numpy as np
def sherpa_model(pr,model_idx,alldists,emissiondelta,inner_radius=False):
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
import pandas as pd
def df2mat(dfdata):
    grids=pd.DataFrame(dfdata.index.str.split('_',expand=True).tolist(), columns=['x','y'], index=dfdata.index)
    grids['x']=pd.to_numeric(grids['x'])
    grids['y']=pd.to_numeric(grids['y'])
    grids['var']=dfdata
    dfmat=grids.pivot(index='y', columns='x', values='var')
    dfmat=dfmat.as_matrix()
    return dfmat

'''
NAME
    Convert a 3D matrix to a multindex dataframe
PURPOSE
    Convert a 3D matrix to a multindex dataframe, where indexes came from a combination of dimnames[0] and dimnames[1] while
    columns come from dimnames[2]
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY

REFERENCES

'''
import numpy as np
import pandas as pd

def array2df(array3D,dimnames):
    #convert it back to Dataframe
    #search which dimensions corresponds to array dimensions
    arrayind=np.argsort(array3D.shape)
    if arrayind[0]==0:
        array3D=np.swapaxes(array3D, 0, 2)
    elif arrayind[1]==0:
        array3D=np.swapaxes(array3D, 1, 2)

    arrayind=np.argsort(array3D.shape)
    if arrayind[0]==1: array3D=np.swapaxes(array3D, 0, 1)

    arrayind=np.argsort(array3D.shape)
    if(arrayind[0]==0 and arrayind[1]==1 and arrayind[2]==2):
     dimnnames_sortednames=[k for k in sorted(dimnames, key=lambda k: len(dimnames[k]))]
     dimnames_len=[len(emissions_array['dimnames'][d]) for  d in dimnnames_sortednames]
     array_new= np.reshape(array3D,(dimnames_len[0]*dimnames_len[1],dimnames_len[2]))
     iterables = [dimnames[dimnnames_sortednames[0]],dimnames[dimnnames_sortednames[1]]]
     df_index=pd.MultiIndex.from_product(iterables, names=[dimnnames_sortednames[0],dimnnames_sortednames[1]])
     dc_df=pd.DataFrame(data=array_new, index=df_index,columns=dimnames[dimnnames_sortednames[2]])
    else:
        print(arrayind)
        exit('something wrong in matrix dimension')
    return dc_df

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
    Denise Pernigotti June 2017
REVISION HISTORY
    From last version of module7_custom.py

REFERENCES

'''
def module7(emissions_nc,concentration_nc, model_nc, intersect_dir,targets_txt,outdir,aggr_zones,outfig='png',rect_coord=None, *progresslog):
    import pandas as pd
    import sys
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from scipy.spatial import distance
    if progresslog:
        progress_dict = read_progress_log(progresslog[0])
        write_netcdf_output = False
    else:
        progress_dict = {'start': 0.0, 'divisor': 1.0}
        write_netcdf_output = True

    if not os.path.exists(outdir):  # added by EPI
        os.makedirs(outdir)

    ################default user definition
    aggr_src=True #set to True for aggregatin sources as in erep
    include_natural=False #include or not natiral PM concentration
    print_areainfo=False #if true check and print which grid points in fua/gcity are actually not modelled
    #############



    grid_txt=intersect_dir+'/grid_intersect'
    gcities_txt=intersect_dir+'/grid_int_gcities'
    fua_txt=intersect_dir+'/grid_int_fua'
    rect_txt=intersect_dir+'/gridint_rect'
    #list of true names for areas IDs
    codes_names=pd.Series({'NUTS_Lv0':'nuts0','NUTS_Lv1':'nuts1','NUTS_Lv2':'nuts2','NUTS_Lv3':'nuts3','FUA_CODE':'fua','GCITY_CODE':'gcities'})
    codes_txt={k: intersect_dir+'/' + codes_names[k] +'_names.txt'for k in codes_names.keys()}
    #list of points used as receptors. Only 'id','lon','lat' and 'Station name' are required columns
    ###############################################

    #read netcdf files, put the data in multindex dataframes and check consistency
    emissions = read_nc(emissions_nc)

    concentration = read_nc(concentration_nc)
    model= read_nc(model_nc)

    if include_natural:
        pmsize=pollutant[-2:]
        salt_nc ='input/pDUST-pSALT',pmsize,'basecase.nc'
        dust_nc = 'input/pDUST-pSALT',pmsize,'basecase.nc'
        salt= read_nc(salt_nc)
        dust= read_nc(dust_nc)

    #check consistency of model and emissions
    if model.loc['coord'].astype(np.float32).equals(emissions.loc['coord'].astype(np.float32)):
        print ('OK latitude and longitude in matrices emissions and model are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")

    #check consistency of model and concentrations
    if model.loc['coord'].astype(np.float32).equals(concentration.loc['coord'].astype(np.float32)):
    #if model.loc['coord'].equals(concentration.loc['coord']):
        print ('OK latitude and longitude in matrices model and conc are the same')
    else:
        sys.exit("latitude and/or longitude are different in loaded matrices")

    #get inner radius and check if the flat weight option is activated
    inner_radius = int(getattr(Dataset(model_nc,'r'), 'Radius of influence'))
    if 'flatWeight' in model.index.levels[0]:
        print('Flatweight approximation is ON with Radius of influence of '+str(inner_radius) +' grid cells')
        if inner_radius >=200:
            sys.exit("Something wrong, radius of influence is >=200 grid cells")
    else:
        if inner_radius is not 200:
            sys.exit("Something wrong, flatweight approximation is deactivated but inner_radius is not 200 grid cells")
        print('Flatweight approximation is OFF')
        inner_radius=False

    #save coord info
    coordinates=emissions.loc['coord',].transpose()
    lat_vec=np.unique(coordinates['lat'])
    lon_vec=np.unique(coordinates['lon'])
    #find grid resolution in degrees
    dlon_res=min(lon_vec[1:]-lon_vec[:-1])
    dlat_res=min(lat_vec[1:]-lat_vec[:-1])

    #grab nuts/areas info from txt
    nuts_info=read_nuts_area(grid_txt,calcall=True)
    nuts_info.update(read_nuts_area(gcities_txt,nullnut='LAND000'))
    nuts_info.update(read_nuts_area(fua_txt,nullnut='LAND000'))
    if(aggr_zones=='rect'):
        nuts_info.update(read_nuts_area('rect',nutsall=nuts_info['ALL_NUTS_Lv0'].copy()))


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

    #grab true names for each area
    area_names_long={k: pd.read_csv(codes_txt[k],delimiter="\t",encoding = "ISO-8859-1",skiprows=1,index_col=0,header=None) for k in codes_names.keys()}
    #reduce string length if needed
    area_names={}
    countries=np.unique(nuts_info['NUTS_Lv0'].index.get_level_values(0))
    for k in area_names_long.keys():
        area_names[k]=area_names_long[k][1].apply(name_short)
    area_names['ALL_NUTS_Lv0']=pd.Series({'ALL_NUTS_Lv0' : 'Europe', 'other' : 'other'},name='EU')
#    area_names['ALL_NUTS_Lv0']=pd.Series({'ALL_NUTS_Lv0' : 'Europe'+str(len(countries)), 'other' : 'other'},name='EU')
    if(aggr_zones=='rect'):
        area_names['rect']=pd.Series({'rect' : 'Domain', 'other' : 'other'},name='EU')
#    totalname='Total'+str(len(countries))
    totalname='Total'
    ########
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


    ###########
    #define the aggregation and increments calculation type depending on aggr_zones
    if aggr_zones=='city':
#        wantedorder=pd.Series( {'3' : 'GCITY_CODE', '2' : 'FUA_CODE', '1' :'NUTS_Lv0','0':'ALL_NUTS_Lv0'})
        wantedorder=pd.Series( {'3' : 'ALL_NUTS_Lv0', '2' : 'NUTS_Lv0', '1' :'FUA_CODE','0':'GCITY_CODE'})
#        wantedorder=pd.Series( {'2' : 'ALL_NUTS_Lv0', '1' : 'NUTS_Lv0', '0' :'FUA_CODE','0':'GCITY_CODE'})
# COM start
    elif aggr_zones=='fuaonly':
        wantedorder=pd.Series( {'1' : 'FUA_CODE', '2':'NUTS_Lv0', '3':'ALL_NUTS_Lv0'})
# COM end
    elif aggr_zones=='nuts':
        wantedorder=pd.Series( {'4' : 'NUTS_Lv3','3' : 'NUTS_Lv2', '2' : 'NUTS_Lv1', '1' :'NUTS_Lv0','0':'ALL_NUTS_Lv0'})
    elif aggr_zones=='rect':
        wantedorder=pd.Series( { '1' :'rect','0':'ALL_NUTS_Lv0'})

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
    dc_ppm={}
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
        #create an emissions delta with zeroes out of the area and emissions detuced by emissions_delta multiplies for the realtive grid area inside the area

        #Core of the sharpa calculations: calculate dc due to each gridpoint for each precursor and macrosector
        start =time.perf_counter()
        dc[idx]=pd.concat(list(map(lambda p: sherpa_model(p, model[idx],dists_array[ix,:],emissions,inner_radius),precursors)))
        dc[idx].index=emissions.index
        end = time.perf_counter()
        #print ("time elapsed ",end-start)

        #aggregate emissions per precursor and area
        nareasurf=pd.concat(list(map(lambda areaname:nuts_info[areaname]['area'].loc[target_info.loc[areaname,'areaid']],target_info.index)),axis=1)
        nareasurf.columns=target_info.index
        nareasurf=nareasurf.fillna(0)
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

         #rename fue elements (same names as city core)
        if 'FUA_CODE' in wantedorder_present['areaid'].values:
            wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']='Commuting zone\n'+wantedorder_present.loc[wantedorder_present['areaid'].values=='FUA_CODE','areaname']
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
            natural=pd.DataFrame(0, index=['Salt','Dust'], columns=dc_inc.columns)
            natural.loc['Dust',totalname]=dust.loc['pDUST-'+pmsize,idx].values*100./concentration[idx]
            natural.loc['Salt',totalname]=salt.loc['pSALT-'+pmsize,idx].values*100./concentration[idx]
            dc_inc=dc_inc.append(natural)
        #plots
        fig={}
        fig[1]=plot_bar(dc_inc,wantedorder_present,totalname)
        plt.close('all')
#        fig[2]=plot_bar(dc_inc_p,wantedorder_present,totalname)
        plt.close('all')
        if emi_units[0] == 'Mg/km2':  # EMA - correct emissions units in plot
            emi_unitsplot = 'ton/year'
        if len(smallareas)>0:
            #fig[3]=plot_dict(dc_inc_flat,idx,coordinates.loc[emissions.columns,])
#            fig[3]=plot_dict(narea_id,idx,coordinates.loc[emissions.columns,])
            plt.close('all')
            for ip,p in enumerate(precursors):
                xlab = ''.join([p, ' emitted mass in ', emi_unitsplot])
                fig[4+ip] = plot_spyder(emi_sum, p)
#                fig[4+ip] = plot_bar(emi_inc.loc[p][smallareas],wantedorder_present.loc[wantedorder_present['areaid'].isin(smallareas)],totalname,plot_opt='noperc',x_label=xlab,leg_loc='upper left')
        plt.close('all')
#        print(dc_inc)
        dc_inc.columns=wantedorder_present['areaname']

        for ids in list(receptors[receptors['target_idx']==idx].index):
            fig[1].savefig(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_sec_bars.'+outfig)
#            fig[2].savefig(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_prec_bars.'+outfig)
            if len(smallareas)>0:
#                fig[3].savefig(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_sec_map.'+outfig)
                for ip,p in enumerate(precursors):
                   if fig[4+ip]:  # added by EMA
                         fig[4+ip].savefig((outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_'+p+'_emi_spys.'+outfig), dpi=300, bbox_inches='tight')
#                    fig[4+ip].savefig(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_'+p+'_emi_bars.'+outfig)
              #dc_inc.to_html(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_total_table.html',classes='table')
            dc_inc.to_csv(outdir+'\\'+ids+'_'+pollutant+'_'+aggr_zones+'_total_table.csv')
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
    #summary.to_html(outdir+'\\AAA_summary_info.html',classes='table')
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
     # pdf from htl template as in http://stackoverflow.com/questions/27387923/combine-existing-figures-into-one-pdf-of-figure-python
    #path_wkthmltopdf = r'C:\Program Files\wkhtmltopdf\bin\wkhtmltopdf.exe'
    #config = pdfkit.configuration(wkhtmltopdf=path_wkthmltopdf)
    #pdfkit.from_file('output/template.html', 'output/output.pdf',configuration=config)
    return emi_sum, wantedorder_present, totalname, narea, dc_inc

if __name__ == '__main__':
    import sys
    import os
    import time
    #from sherpa_auxiliaries import write_progress_log

    # module 7  test inputs
    module = 7
    # if it doesn't exist strart=0 and dividsor=1
    progresslog = 'input/progress.log'

    # run module 1 without progress log
    ############################################### user input data
#    sherpa_version='inputFlatWeightChimere_7km_nuts'
    sherpa_version='20170322_v18_SrrResults_PotencyBased'
    testarea='Pariscity' #may be any area as long as the file testarea_targets.txt is present in input, contains a list of lat/lon
    pollutant='PM25' #may be 'PM25' or 'PM10' or NOx
    aggr_zones='city' #may be 'city','nuts' or 'rect' (in this case the domain defined with ll and ur)
    #rect_coord={'ll':{'lat':47.9375,'lon':-2.2500},'ur':{'lat':53.0000,'lon':6.3750}}
    outfig='jpg' #'pnd' of 'pdf'
    ###############################################
    ############################################### input files
    pollconc=pollutant+'_Y'
    pollmodel=pollutant+'_Y'
    if pollutant is 'NOx':
        pollutant='NO2'
        pollconc='NO2_NO2eq_Y_mgm3'
        pollmodel='NO2eq_Y'
    elif pollutant not in ['PM10','PM25']:
        sys.exit(pollutant +'is not implemented in this module')
    emissions='input/'+sherpa_version+'/1_base_emissions/BC_emi_'+pollutant+'_Y.nc'
    concentration='input/'+sherpa_version+'/2_base_concentrations/BC_conc_'+pollconc+'.nc'
    model='input/'+sherpa_version+'/3_source_receptors/SR_'+pollmodel+'_20170322_potencyBased.nc'   # changed ema** #info on areas and percentage of grids in areas
    selection_dir='input/selection'
    target_list='input/'+testarea+'_targets.txt'
    outdir='output/'+sherpa_version+'/'+testarea
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #directory for emissions and SR relationship only
    print('SHERPA input will be searched in '+model)
    print('with emissions in '+emissions)
    print('with BC concentrations in conc variable in '+concentration)
    print('while txt files for grid intersect nuts, fua etc will be searced in local input directory '+selection_dir)

    # run module 7 with progress log
    #proglog_filename = 'output/proglog'
    start = time.perf_counter()
#    out_dc=module7(emissions,concentration, model,selection_dir,target_list,outdir,aggr_zones,outfig)
    [emi_sum, wantedorder_present, totalname, narea, dc_inc] = module7(emissions,concentration, model,selection_dir,target_list,outdir,aggr_zones,outfig)
    stop = time.perf_counter()
    print('Module 7 run time: %s sec.' % (stop-start))
    #remove(proglog_filename)

    pass