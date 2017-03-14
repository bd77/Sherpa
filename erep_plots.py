import numpy as np
import pandas as pd
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
def plot_dict(dc_dic,idx,lats,lons,grid_grads=3,unit='%'):  
    lat_0=lats.loc[idx]
    lon_0=lons.loc[idx]
    ll_lon=lon_0-grid_grads
    ll_lat=lat_0-grid_grads
    ur_lon=lon_0+grid_grads
    ur_lat=lat_0+grid_grads
    lon, lat = np.meshgrid(lons.unique(),lats.unique())    
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
    clevs = [5,10,15,20,25,30,35,40,45,50,55,60]
    cs = m.contourf(x,y,dfmat,clevs);
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label(unit)
    #drow target
    x,y = m(lon_0, lat_0)
    m.plot(x, y, 'x', markersize=6,color='k',fillstyle='none')
    #plt.title(idx)
    #plt.savefig(fileout)
    #plt.close(fig)
    #plt.show()  
    return fig

'''
NAME
   Stacked barplot using plt
PURPOSE
    Stacked barplot using plt, http://chrisalbon.com/python/matplotlib_stacked_bar_plot.html
    Input parameter ftx is the font size of labels, 16 is the default
PROGRAMMER(S)
    Denise Pernigotti 
REVISION HISTORY
    
REFERENCES
    
'''
def plot_bar(dfdata,varplot,ftx=16):  
    dfdata=dfdata[varplot['areaid']].transpose()
    varnames=dfdata.columns
    areanames=varplot['areaname']
    areanames.index=varplot['areaid']
    varplot_left=pd.Series(np.zeros(len(dfdata.index)),index=dfdata.index)
    addsum=dfdata.sum(axis=1)
    for ind in range(1,len(dfdata.index)):
        varplot_left[ind]=varplot_left[ind-1]+addsum[dfdata.index[ind-1]]
    varplot_left['total']=0.
    colors={'PPM':'grey','SOx':'#ff0000','NOx': '#0000ff','NH3':'#b22222','NMVOC':'#8470ff',
    'Traffic': '#0000ff','Industry':'#ff0000','Agriculture':'#b22222','Residential':'#32cd32',
    'Offroad':'#8470ff','Other':'#bebebe','Salt':'#ccffe5','Dust':'#ffffcc'}
    # Create the general blog and the "subplots" i.e. the bars
    plt.close('all')
    f, ax1 = plt.subplots(1, figsize=(10,5))
    ax1.set_xlim([0, 100])
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(ftx) 
    # Set the bar width
    bar_width = 0.75

    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(varplot))] 

    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l] 

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
        #alpha=0.5, 
        # with color
        color=colors[varnames[0]])

    # Create a bar plot, in position bar_1
    for ivar in range(1,len(varnames)):
        ax1.barh(bar_l, 
        # using the mid_score data
        dfdata[varnames[ivar]], 
        # set the width
        height=bar_width,
        # with pre_score and mid_score on the bottom
        left=list(varplot_left+dfdata[varnames[range(0,ivar)]].sum(axis=1)), 
        # with the label post score
        label=varnames[ivar], 
        # with alpha 0.5
        #alpha=0.5, 
        # with color
        color=colors[varnames[ivar]])

# set the y ticks with names
    plt.yticks(tick_pos, areanames,size=12)

# Set the label and legends
    ax1.set_xlabel("% of total mass",fontsize=ftx)
    #ax1.set_ylabel("Areas")
    plt.legend(loc='lower right')

# Set a buffer around the edge
    plt.ylim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
    return f


