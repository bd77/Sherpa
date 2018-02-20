# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:15:49 2017

@author: peduzem
"""
import numpy as np
import matplotlib as mpl
import matplotlib.patches as mpatches

def figsize(scale):
    fig_width_pt = 390                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt
#def savefig(filename):
#    plt.savefig('{}.pgf'.format(filename))
#    plt.savefig('{}.pdf'.format(filename))
#, bbox_inches = "tight", dpi=300


# for importing matlab files
import scipy.io as sio
from PIL import Image
from osgeo import gdal, ogr, osr  #conda install -c conda-forge gdal
from netCDF4 import  Dataset # for using netCDF files
from time import time  # for module1
import os as os
from os import remove
import pandas as pd  # conda install pandas


#from matplotlib import rc
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)
#import plotly.plotly as py
import matplotlib.lines as mlines
import matplotlib.markers as mmarks
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter,ScalarFormatter
import matplotlib.patches as mpatches
# for plotting
from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap

from healthia import calc_impacts
from sherpa_auxiliaries import (create_emission_reduction_dict,
                                create_emission_dict, create_window,
                                read_progress_log, write_progress_log)
from sherpa_auxiliaries_epe import (gridint_toarray, write_nc, tiftogridgeneral,
                                    save_obj, load_obj, bm_precompute)
from sherpa_globals import (path_base_conc_cdf_test, path_salt_conc_cdf_test,
                            path_dust_conc_cdf_test, path_pop_mat_test,
                            path_mortbaseline, path_emission_cdf_test,
                            path_model_cdf_test, path_result_cdf_test)
from module7_custom_ep import read_nuts_area

#def plot_secregeffect(precursor_lst, name, sources, path_results, path_figures):
#    """ Plots aggregated sectors efficiency ratios in FUA vs the efficiency
#    ratio of only the fua
#    """
#    if not os.path.exists(path_figures):
#                os.makedirs(path_figures)
#    N = len(sources)
#    ind = np.arange(N)
#    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
#                 'other']
#    # for the plots, marks corresponding to precursors
#    marks = ['o', 'v', 's', '*']
#    # Plots by sectors and FUA
#
#
#    fig, ax = plt.subplots(figsize=(14, 8))
#    ax.minorticks_on()
#    patches = []
#    for ind, precursor in enumerate(precursor_lst):
#        df_res = pd.read_csv(path_results +
#                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
#        # Dataframe with reduction per sectors and per FUA
#        df_sr = pd.read_csv(path_results +
#                        'df_sr_{}.csv'.format(precursor), index_col=[0, 1])
#        df_scaty = pd.DataFrame(index=sources, columns=sect_aggr)
#        df_scatx = pd.DataFrame(index=sources, columns=(['ref_ratio', 'fua_ratio', 'fua_ref']))
#        df_scatx['ref_ratio'] = (df_res[name].ix[:, 'country'] /
#                    (df_res['emi_reg'].ix[:, 'country']))
#        df_scatx['fua_ratio'] = (df_res[name].ix[:, 'FUA'] /
#                    (df_res['emi_reg'].ix[:, 'FUA']))
#        df_scatx['fua_ref'] = df_scatx['fua_ratio'] / df_scatx['ref_ratio']
#        for inds, sec in enumerate(sect_aggr):
#            df_scaty[sec] = ((df_sr[name].loc[:, sec] /
#                             (df_sr['emi_sec'].loc[:, sec])) /
#                             df_scatx['ref_ratio'])
#            plt.plot(df_scatx['fua_ref'], df_scaty[sec],
#                     marks[ind], label = '{}'.format(sec), color= 'C{}'.format(inds))
#
#    for ind, sec in enumerate(sect_aggr):
#        patches.append(mpatches.Patch(color= 'C{}'.format(ind), label='{}'.format(sec)))
#
#    for ind, prec in enumerate(precursor_lst):
#        marker = mlines.Line2D([], [], color='black', marker=marks[ind],
#                          markersize=6, label=precursor_lst[ind], linestyle = 'None')
#        patches.append(marker)
#
#    plt.xlabel('FUA eff. ratio, $\\rho_p^{{FUA}}$, for {}'.format(name))
#    plt.ylabel('FUA sec. eff. ratio, $\\rho_p^{{FUA, MS}}$, for {}'.format(name))
#    plt.legend(handles=patches)
#
#    fig.savefig(path_figures + 'perprecpersecFUA_{}.png'.format(name), bbox_inches = "tight", dpi=300)


def plot_effpot(precursor_lst, name, sources, path_results,
                  path_figures):
    # aggregated secotrial and per precursor efficiency and potential
#    if yvar == 'effratio':
#        yaxiskword = 'Efficiency ratio, $\\rho_p^{{MS}}$'
#    yvar = 'eff'

    ### Part 1 - sectors contribution:
    yaxiskword = 'Potency, $\\eta_{{p,sector}}$'
    marks = ['o', 'v', 's', '*']
    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
                 'other']
    # Plots by sectors:
    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()

    # base case dataframe (concentration, exposure, emi)
    df_bc = pd.read_csv(path_results +
                        'df_bc.csv', index_col=[0])
    patches = []
    for indp, precursor in enumerate(precursor_lst):
        # result of regional reductions (FUA, not_FUA, country)
#        df_res = pd.read_csv(path_results +
#                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
#        df_sr = pd.read_csv(path_results +
#                        'df_sr_{}.csv'.format(precursor), index_col=[0, 1])
        # result of sectorial reductions by country (concentration change give precred)
        df_sec = pd.read_csv(path_results +
                        'df_sec_{}.csv'.format(precursor), index_col=[0, 1])
#        df_plots =  pd.read_csv(path_results +
#                        'df_plots_{}{}.csv'.format(name, precursor),
#                                      index_col=[0])
#        # dataframe with x axis : potential
        df_scatx = pd.DataFrame(index=sources, columns=sect_aggr)
        # dataframe with y axix : efficiency (potency)
        df_scaty = pd.DataFrame(index=sources, columns=sect_aggr)
        aggsecemi = pd.DataFrame(index=sources, columns=sect_aggr)
        marker = mlines.Line2D([], [], color='black', marker=marks[indp],
                          markersize=6, label=precursor_lst[indp], linestyle = 'None')
        patches.append(marker)
        for ind, sec in enumerate(sect_aggr):
            if sec == 'industry':
                # concentration change per precred
                aggsec = (df_sec[name].loc[:, 'MS1'] +
                          df_sec[name].loc[:, 'MS3'] +
                          df_sec[name].loc[:, 'MS4'])
                # emission of the sector
                aggsecemi[sec] = (df_sec['emi_sec'].loc[:, 'MS1'] +
                          df_sec['emi_sec'].loc[:, 'MS3'] +
                          df_sec['emi_sec'].loc[:, 'MS4'])
                df_scatx[sec] = ((aggsec / df_bc[name]))
                df_scaty[sec] = (aggsec / (aggsecemi[sec]))
            elif sec == 'residential':
                df_scatx[sec] = (((df_sec[name].loc[:, 'MS2']) / df_bc[name]))
                aggsecemi[sec] =df_sec['emi_sec'].loc[:, 'MS2']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS2']) / (df_sec['emi_sec'].loc[:, 'MS2'])
            elif sec == 'agriculture':
                df_scatx[sec] = (((df_sec[name].loc[:, 'MS10']) / df_bc[name]))
                aggsecemi[sec] = df_sec['emi_sec'].loc[:, 'MS10']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS10']) / (df_sec['emi_sec'].loc[:, 'MS10'])
            elif sec == 'transport':
                df_scatx[sec] = (((df_sec[name].loc[:, 'MS7']) / df_bc[name]))
                aggsecemi[sec]=df_sec['emi_sec'].loc[:, 'MS7']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS7']) / (df_sec['emi_sec'].loc[:, 'MS7'])
            elif sec == 'other':
                aggsec =(df_sec[name].loc[:, 'MS5'] +
                                 df_sec[name].loc[:, 'MS6'] +
                                 df_sec[name].loc[:, 'MS8'] +
                                 df_sec[name].loc[:, 'MS9'])
                aggsecemi[sec] =(df_sec['emi_sec'].loc[:, 'MS5'] +
                                 df_sec['emi_sec'].loc[:, 'MS6'] +
                                 df_sec['emi_sec'].loc[:, 'MS8'] +
                                 df_sec['emi_sec'].loc[:, 'MS9'])
                df_scatx[sec] = ((aggsec / df_bc[name]))
                df_scaty[sec] = (aggsec / (aggsecemi[sec]))

            plt.plot(df_scatx[sec], df_scaty[sec], marks[indp], label = '{}'.format(sec), color= 'C{}'.format(int(ind)))
        df_scatx.to_csv(path_figures + 'df_potential_aggr_{}{}.csv'.format(name, precursor))
        df_scaty.to_csv(path_figures + 'df_efficiency_aggr_{}{}.csv'.format(name, precursor))
        aggsecemi.to_csv(path_figures + 'df_aggsecemi_aggr_{}.csv'.format(precursor))
            #plt
    plt.xlabel('potential, $\phi_p^{{Sector}}$, for {}'.format(name))
    plt.ylabel('{}, for {}'.format(yaxiskword, name))

    for ind, sec in enumerate(sect_aggr):
        patches.append(mpatches.Patch(color= 'C{}'.format(int(ind)), label='{}'.format(sec)))
    plt.legend(handles=patches)
    plt.show()
    fig.savefig(path_figures + 'perprecpersectex_{}.pgf'.format(name)) #,
    fig.savefig(path_figures + 'perprecpersectex_{}.pdf'.format(name))


    ### Part 2 - areas contribution:
    yaxiskword = 'potency, $\\eta_{{p,area}}$'
    marks = ['o', 'v', 's', '*']
    areas = ['FUA', 'not_FUA']
    # Plots by sectors:
    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()

    # base case dataframe (concentration, exposure, emi)
    df_bc = pd.read_csv(path_results +
                        'df_bc.csv', index_col=[0])
    patches = []
    for indp, precursor in enumerate(precursor_lst):
        marker = mlines.Line2D([], [], color='black', marker=marks[indp],
                          markersize=6, label=precursor_lst[indp], linestyle = 'None')
        patches.append(marker)
        # result of regional reductions (fua, not_fua, country)
        df_scatx = pd.DataFrame(index=sources, columns=areas)
        # dataframe with y axix : efficiency (potency)
        df_scaty = pd.DataFrame(index=sources, columns=areas)
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_scaty['FUA'] = (df_res[name].ix[:, 'FUA'] /
             df_res['emi_reg'].ix[:, 'FUA'])

        df_scaty['not_FUA'] = (df_res[name].ix[:, 'not_FUA'] /
            (df_res['emi_reg'].ix[:, 'not_FUA']))

        df_scatx['FUA'] = ((df_res[name].ix[:, 'FUA'] / df_bc[name]))
        df_scatx['not_FUA'] = ((df_res[name].ix[:, 'not_FUA'] / df_bc[name]))
        df_scatx.to_csv(path_figures + 'df_potential_area_{}{}.csv'.format(name, precursor))
        df_scaty.to_csv(path_figures + 'df_efficiency_area_{}{}.csv'.format(name, precursor))

        plt.plot(df_scatx['FUA'], df_scaty['FUA'], marks[indp], color= 'C{}'.format(int(1)))
        plt.plot(df_scatx['not_FUA'], df_scaty['not_FUA'], marks[indp], color= 'C{}'.format(int(2)))
    plt.xlabel('potential, $\phi_{{p,area}}$, for {}'.format(name))
    plt.ylabel('{}, for {}'.format(yaxiskword, name))
    patches.append(mpatches.Patch(color= 'C{}'.format(int(1)), label = '{}'.format('FUA')))
    patches.append(mpatches.Patch(color= 'C{}'.format(int(2)), label = '{}'.format('non FUA')))
    plt.legend(handles=patches)
    plt.show()
    fig.savefig(path_figures + 'perprecperareatex_{}.pgf'.format(name)) #,
    fig.savefig(path_figures + 'perprecperareatex_{}.pdf'.format(name))


def plot_ratios_fuas(name, sources, precursor_lst, path_results, path_figures):
#    marks = ['-o', '-v', '-s', '-*']
    marks = ['o', 'v', 's', '*']
    if not os.path.exists(path_figures):
                os.makedirs(path_figures)
    targets = ['ref_ratio', 'fua_ratio', 'not_fua_ratio']
    df_plots = pd.DataFrame(index=sources, columns=targets)
    for indp, precursor in enumerate(precursor_lst):
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_plots['ref_ratio'] = (df_res[name].ix[:, 'country'] /
            (df_res['emi_reg'].ix[:, 'country']))

        df_plots['fua_ratio'] = (df_res[name].ix[:, 'FUA'] /
            (df_res['emi_reg'].ix[:, 'FUA']))

        df_plots['not_fua_ratio'] = (df_res[name].ix[:, 'not_FUA'] /
            (df_res['emi_reg'].ix[:, 'not_FUA']))

        N = len(sources)
        ind = np.arange(N)
        plt.clf()
        fig= plt.figure(figsize=figsize(1))
        ax = fig.add_subplot(111)
        ax.minorticks_on()
#        fig, ax = plt.subplots(figsize=(14, 8))
        ax.set_ylim([0,3])
#        ax = plt.minorticks_on()
        ax.tick_params(axis='x',which='minor',bottom='off')
        plt.xticks(ind, sources)
        p1 = plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'], '-', color = 'black')
        p2 = plt.plot(ind, df_plots['fua_ratio']/df_plots['ref_ratio'], marks[indp])
        p3 = plt.plot(ind, df_plots['not_fua_ratio']/df_plots['ref_ratio'], marks[indp])
        plt.ylabel('{} performance ratio, $\\rho_{{s,{}}}$'.format(name, precursor))
        plt.legend((p1[0], p2[0], p3[0]), ('uni. red.', 'FUA', 'Non FUA'), loc='best', ncol=1)
#        plt.savefig(path_figures + 'regionalratios{}_{}.png'.format(name, precursor),
#                    bbox_inches = "tight", dpi=300)
        fig.savefig(path_figures + 'regionalratiostex{}_{}.pgf'.format(name, precursor))
        fig.savefig(path_figures + 'regionalratiostex{}_{}.pdf'.format(name, precursor))
        plt.show()

if __name__ == '__main__':

    # -------------------------------------------------------------------------
    # Display results
    # -------------------------------------------------------------------------


    sources = ['BE', 'DE','ES', 'FR', 'IT', 'UK']
#    path_results = 'spaflextest\\'

    path_results = 'spaflextest3\\'

    path_figures = path_results
    name = 'concentration' # have not checked if everything is correct for concentration
    precursor_lst = ['PPM', 'SOx', 'NOx', 'NH3']#, 'NMVOC']#, 'NMVOC']
#    rootgrp = Dataset(path_model_cdf_test, 'r')
#    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
#    rootgrp.close()

    path_surf_nc = 'input/JRC01.nc'
    # Area of each cell in the domain
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    path_tiff = 'input/pop/7km_Qi_2010.tif'
    popall = tiftogridgeneral(path_tiff)
    level = 'NUTS_Lv0'
    parea = 'parea'
    colfuas = ['perpop','perarea','peremi'] #percentage population , percentage area
    df_fuas = pd.DataFrame(index=sources, columns=colfuas)
    precursor = 'PPM'
    df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
    for country in sources:
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
#        if not os.path.exists(path_areacountry_nc):
#            areacountry = gridint_toarray('NUTS_Lv0', 'parea', country)
#            write_nc(areacountry, path_areacountry_nc, 'AREA', '%')
#        else:
        print(('WARNING: Loading previously calculated country areas'))
        rootgrp = Dataset(path_areacountry_nc, mode='r')
        areacountry = rootgrp.variables['AREA'][:]
            # bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()
        path_areafuas_nc = ('asstest\\workdirreg\\area_fuas_{}.nc'.format(country))
        rootgrp = Dataset(path_areafuas_nc, mode='r')
        areafuas = rootgrp.variables['AREA'][:]
        df_fuas.loc[country]['perpop']= np.array(popall*areafuas).sum()/np.array(popall*areacountry).sum()*100
        df_fuas.loc[country]['perarea']= np.array(surf*areafuas).sum()/np.array(surf*areacountry).sum()*100
        df_fuas.loc[country]['peremi']= df_res['emi_reg'].ix[country, 'FUA']/df_res['emi_reg'].ix[country, 'country']*100
        # aggregated secotrial and regional effect
#    plot_secregeffect(precursor_lst, name, sources, path_results,
#                  path_figures)
    plot_ratios_fuas(name, sources, precursor_lst, path_results, path_figures)

    # aggregated secotrial and efficiency and potential
    yvar= 'eff' # or 'effratio'
    plot_effpot(precursor_lst, name, sources, path_results,
                path_figures)

    # ms_list = ['MS1', 'MS2', 'MS7', 'MS10']
    # PLOT RATIOS for aggregated sectors
    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
                 'other']
    targets = ['ref_ratio']
    N = len(sources)
    ind = np.arange(N)
    marks = ['o', 'v', 's', '*']
#    marks = ['-o', '-v', '-s', '-*']
    for indp, precursor in enumerate(precursor_lst):
        df_res = pd.read_csv(path_results + 'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_plots = pd.DataFrame(index=sources, columns=targets)
        df_plots['ref_ratio'] = (df_res[name].ix[:, 'country'] /
             ((df_res['emi_reg'].ix[:, 'country'])))
        df_sec = pd.read_csv(path_results + 'df_sec_{}.csv'.format(precursor), index_col=[0, 1])
        for sec in sect_aggr:
            if sec == 'industry':
                # concentration change per precred
                aggsec = (df_sec[name].loc[:, 'MS1'] +
                          df_sec[name].loc[:, 'MS3'] +
                          df_sec[name].loc[:, 'MS4'])
                # emission of the sector
                aggsecemi = (df_sec['emi_sec'].loc[:, 'MS1'] +
                          df_sec['emi_sec'].loc[:, 'MS3'] +
                          df_sec['emi_sec'].loc[:, 'MS4'])
                df_plots[sec] = (aggsec / (aggsecemi))
            elif sec == 'residential':
                df_plots[sec] = (df_sec[name].loc[:, 'MS2']) / (
                                 df_sec['emi_sec'].loc[:, 'MS2'])
            elif sec == 'agriculture':
                df_plots[sec] = (df_sec[name].loc[:, 'MS10']) / (
                                 df_sec['emi_sec'].loc[:, 'MS10'])
            elif sec == 'transport':
                df_plots[sec] = (df_sec[name].loc[:, 'MS7']) / (
                                 df_sec['emi_sec'].loc[:, 'MS7'])
            elif sec == 'other':
                aggsec =(df_sec[name].loc[:, 'MS5'] +
                                 df_sec[name].loc[:, 'MS6'] +
                                 df_sec[name].loc[:, 'MS8'] +
                                 df_sec[name].loc[:, 'MS9'])
                aggsecemi =(df_sec['emi_sec'].loc[:, 'MS5'] +
                                 df_sec['emi_sec'].loc[:, 'MS6'] +
                                 df_sec['emi_sec'].loc[:, 'MS8'] +
                                 df_sec['emi_sec'].loc[:, 'MS9'])
                df_plots[sec] = (aggsec / (aggsecemi))

        df_plots.to_csv(path_figures + 'df_plots_aggr_{}{}.csv'.format(name, precursor))
        plt.clf()
        fig= plt.figure(figsize=figsize(1))
        ax = fig.add_subplot(111)
        ax.minorticks_on()
#        fig, ax = plt.subplots(figsize=(14, 8))
#        ax.set_yscale('log', basey=2)
        ax.set_ylim([0,3])
#        ax.set_yticks(np.arange(0,4))
#        ax.get_yaxis().set_major_formatter(ScalarFormatter())

#        ax = plt.minorticks_on()
        ax.tick_params(axis='x',which='minor',bottom='off')
        plt.xticks(ind, sources)
        plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'],'-' , color = '#000000', label = 'uni. red.') #marks[indp]
        for inds, sec in enumerate(sect_aggr):
            plt.plot(ind, df_plots['{}'.format(sec)]/df_plots['ref_ratio'],
                 marks[indp], label='{}'.format(sec),
                 color='C{}'.format(inds))

        plt.ylabel('{} performance ratio, $\\rho_{{s,{}}}$'.format(name, precursor))
        plt.legend(ncol=2, loc='best')

        
        fig.savefig(path_figures + 'sectorialratiosnew{}_{}.png'.format(name, precursor),
                    bbox_inches = "tight", dpi=300)
        fig.savefig(path_figures + 'sectorialratiostex{}_{}.pgf'.format(name, precursor))
        fig.savefig(path_figures + 'sectorialratiosex{}_{}.pdf'.format(name, precursor))

        plt.show()

#    color= 'C{}'.format(int(ms[2:])-1)
#    sources = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES',
#               'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
#               'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
#
##    sources = ['EL']
#
#    path_results = 'asstest\\'
#    precred = 50
#    N = len(sources)
#    ind = np.arange(N)
#    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
#    # for name in names:
#    name = 'exposure'
#    precursor_lst = ['PPM']#, 'SOx', 'NOx', 'NH3']#, 'NMVOC']
#    marks = ['o', 'v', 's', '*']
#
#    # Plots by sectors:
#    fig, ax = plt.subplots(figsize=(14, 8))
#    ax = plt.minorticks_on()
#    patches = []
#    for ind, precursor in enumerate(precursor_lst):
#        df_bc = pd.read_csv(path_results +
#                        'df_bc.csv', index_col=[0])
#        df_res = pd.read_csv(path_results +
#                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
#        df_sec = pd.read_csv(path_results +
#                        'df_sec_{}.csv'.format(precursor), index_col=[0, 1])
#        df_plots =  pd.read_csv(path_results +
#                        'df_plots_{}{}.csv'.format(name, precursor),
#                                      index_col=[0])
#        df_scatx = pd.DataFrame(index=sources, columns=ms_list)
#        df_scaty = pd.DataFrame(index=sources, columns=(['ref_ratio', 'fua_ratio', 'fua_ref']))
#        df_scaty['ref_ratio'] = (df_res[name].ix[:, 'country'] /
#                    (df_res['emi_reg'].ix[:, 'country'] * precred / 100))
#        df_scaty['fua_ratio'] = (df_res[name].ix[:, 'FUA'] /
#                    (df_res['emi_reg'].ix[:, 'FUA'] * precred / 100))
#        df_scaty['fua_ref'] = df_scaty['fua_ratio'] / df_scaty['ref_ratio']
#        for ms in ms_list:
#            # df_scaty[ms] = df_sec[name].ix[:, ms]/df_res[name].ix[:, 'country']
##            df_scaty[ms] = (df_sec[name].ix[:, ms] * 1/ (precred / 100) )/df_bc[name]
#            df_scatx[ms] =  df_plots['{}'.format(ms)]/df_plots['ref_ratio']
#            plt.plot(df_scatx[ms], df_scaty['fua_ref'], marks[ind], label = '{}'.format(ms))
#
#    for ind, ms in enumerate(ms_list):
#        patches.append(mpatches.Patch(color= 'C{}'.format(ind), label='{}'.format(ms)))
#
#    for ind, ms in enumerate(precursor_lst):
#        marker = mlines.Line2D([], [], color='black', marker=marks[ind],
#                          markersize=6, label=precursor_lst[ind], linestyle = 'None')
#        patches.append(marker)
#
#    plt.ylabel('FUA ratio, $\\rho_p^{{FUA}}$, for {}'.format(name))
#    plt.xlabel('Efficiency ratio, $\\rho_p^{{MS}}$, for {}'.format(name))
#    plt.legend(handles=patches)
#    plt.savefig(path_results + 'perprecpersecPPM_{}.png'.format(name),
#                    bbox_inches = "tight", dpi=300)
#    # plt.legend(markerscale=0.7, scatterpoints=1)
#    plt.show()


##########
##########
##########




##########
##########
##########

