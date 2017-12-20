# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:15:49 2017
Script used to produce the data for the article
" Impact of source-receptor model spatial flexibility on AQ policies"
@author: peduzem
"""

# for importing matlab files
#import scipy.io as sio
#from PIL import Image
#from osgeo import gdal, ogr, osr  #conda install -c conda-forge gdal
from netCDF4 import  Dataset # for using netCDF files
import numpy as np  # for scientific operators
from time import time  # for module1
from os import remove
import os as os
import pandas as pd  # conda install pandas
# for plotting
# from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap
import matplotlib as mpl
import matplotlib.pyplot as plt
#import plotly.plotly as py
import matplotlib.lines as mlines
#import matplotlib.markers as mmarks
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.patches as mpatches

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
from module7_custom import read_nuts_area

import module1 as shrp

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

def precompute(path_areared_nc):
    """
    Precompute function: computes given an area the correpsonding emission:
    Input:
        - path_areared_nc: path to the area where the reduction is applied
    Output:
        - emi: total emission per pollutant [Mg]
        - emi_sec: total emissions  per pollutant per sector [Mg]
    Files needed:
        - nc of the area of each country in:
            'workdir\\area_{}.nc'.format(country)

    Created on Thu Apr 20 16:26:42 2017
    @author: peduzem
    """
    # Total emissions by country and precursor
    # -------------------------------------------------------------------------
    # create a dictionary with emissions per precursor, macrosector and postion
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()
    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)
    # Area of each cell in the domain
    path_surf_nc = 'input/JRC01.nc'
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    # area_units = rootgrp.variables['surface'].units
    rootgrp.close()

    # Percentage of each cell in the domain
    rootgrp = Dataset(path_areared_nc, mode='r')
    area = rootgrp.variables['AREA'][:]
    # area_units = rootgrp.variables['AREA'].units
    rootgrp.close()

    emi = {precursor: np.sum(
           np.sum((emission_dict[precursor]*surf * area/100), axis=0))
           for precursor in precursor_lst}

    emi_sec = {(precursor, ms):
               np.sum(emission_dict[precursor][ms]*surf * area/100)
               for precursor in precursor_lst for ms in np.arange(0, 10)}

    return emi, emi_sec

def plot_ratios_aggsec(name, prec_lst, sources, path_figures, path_results):
    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
                 'other']
    targets = ['ref_ratio']
    N = len(sources)
    ind = np.arange(N)
    marks = ['-o', '-v', '-s', '-*','-D']
    for indp, precursor in enumerate(prec_lst):
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
                df_plots[sec] = (aggsec.where(aggsec > thr) / (aggsecemi))
            elif sec == 'residential':
                df_plots[sec] = df_sec[name].loc[:, 'MS2'].where(df_sec[name].loc[:, 'MS2']>thr) / (
                                 df_sec['emi_sec'].loc[:, 'MS2'])
            elif sec == 'agriculture':
                df_plots[sec] = (df_sec[name].loc[:, 'MS10'].where(df_sec[name].loc[:, 'MS10']>thr)) / (
                                 df_sec['emi_sec'].loc[:, 'MS10'])
            elif sec == 'transport':
                df_plots[sec] = (df_sec[name].loc[:, 'MS7'].where(df_sec[name].loc[:, 'MS7'] >thr)) / (
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
                df_plots[sec] = aggsec.where(aggsec> thr) / aggsecemi

        df_plots.to_csv(path_figures + 'df_plots_aggr_{}{}.csv'.format(name, precursor))
        plt.clf()
        fig= plt.figure(figsize=figsize(1))
        ax = fig.add_subplot(111)
        ax.minorticks_on()
#        fig, ax = plt.subplots(figsize=(14, 8))
        ax.set_ylim([0,3])
#        ax = plt.minorticks_on()
        ax.tick_params(axis='x',which='minor',bottom='off')
        plt.xticks(ind, sources)
        plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'], marks[indp], color = '#000000', label = 'country ref')
        for inds, sec in enumerate(sect_aggr):
            plt.plot(ind, df_plots['{}'.format(sec)]/df_plots['ref_ratio'],
                 marks[indp], label='{}'.format(sec),
                 color='C{}'.format(inds))

        plt.ylabel('{} performance ratio, $\\rho_{{s,{}}}$'.format(name, precursor))
        plt.legend(ncol=2, loc='best')
        fig.savefig(path_figures + 'sectorialratiosagg{}_{}.png'.format(name, precursor),
                    bbox_inches = "tight", dpi=300)
        fig.savefig(path_figures + 'sectorialratiosagg{}_{}.pgf'.format(name, precursor))
        fig.savefig(path_figures + 'sectorialratiosagg{}_{}.pdf'.format(name, precursor))
        plt.show()

def plot_ratios_secs(name, prec_lst, sources, path_figures, path_results):
    N = len(sources)
    ind = np.arange(N)
    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
    targets = ['ref_ratio']
    marks = ['-o', '-v', '-s', '-*', '-D']
#    prec_lst = ['NMVOC']
#    names = ['exposure']
    for indp, prec in enumerate(prec_lst):
        df_sec = pd.read_csv(path_results +
                             'df_sec_{}.csv'.format(prec), index_col=[0, 1])
        df_res = pd.read_csv(path_results +
                             'df_res_{}.csv'.format(prec), index_col=[0, 1])
        for name in names:
            df_plots = pd.DataFrame(index=sources, columns=targets)
            df_plots['ref_ratio'] = (df_res[name].ix[:, 'country'] /
                     ((df_res['emi_reg'].ix[:, 'country'])))
            for ms in ms_list:
                df_plots['{}'.format(ms)] = (df_sec[name].ix[:, ms].where(df_sec[name].ix[:, ms]>thr) /
                     ((df_sec['emi_sec'].ix[:, ms])))
            df_plots.to_csv(path_results + 'df_plots_{}{}.csv'.format(name, prec))
#            fig, ax = plt.subplots(figsize=(14, 8))
            fig= plt.figure(figsize=figsize(1))
            ax = fig.add_subplot(111)
            ax.set_ylim([0,6])
            ax = plt.minorticks_on()
            ax = plt.tick_params(axis='x',which='minor',bottom='off')
            plt.xticks(ind, sources)
            plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'], marks[indp], color = '#000000', label = 'country ref')
            for ms in ms_list:
                plt.plot(ind, df_plots['{}'.format(ms)]/df_plots['ref_ratio'],
                         marks[indp], label='{}'.format(ms),
                         color='C{}'.format(int(ms[2:])-1))
            plt.ylabel('{} performance ratio, $\\rho_{{s,{}}}$'.format(name, prec))
            plt.legend(ncol=2, loc='best')
            plt.savefig(path_results + 'sectorialratios{}_{}.png'.format(name, prec),
                        bbox_inches = "tight", dpi=300)
            fig.savefig(path_figures + 'sectorialratios{}_{}.pgf'.format(name, prec))
            plt.show()



def plot_effpot(prec_lst, name, sources, path_results,
                  path_figures):
    # aggregated secotrial and per precursor efficiency and potential
#    if yvar == 'effratio':
#        yaxiskword = 'Efficiency ratio, $\\rho_p^{{MS}}$'
#    yvar = 'eff'

    ### Part 1 - sectors contribution:
    yaxiskword = 'Efficiency, $\\eta_{{p,sector}}$'
    marks = ['o', 'v', 's', '*', 'D']
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
    for indp, precursor in enumerate(prec_lst):
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
                          markersize=6, label=prec_lst[indp], linestyle = 'None')
        patches.append(marker)
        for ind, sec in enumerate(sect_aggr):
            if sec == 'industry':
                # concentration contribution
                aggsec = (df_sec[name].loc[:, 'MS1'] +
                          df_sec[name].loc[:, 'MS3'] +
                          df_sec[name].loc[:, 'MS4'])
                # emission of the sector
                aggsecemi[sec] = (df_sec['emi_sec'].loc[:, 'MS1'] +
                          df_sec['emi_sec'].loc[:, 'MS3'] +
                          df_sec['emi_sec'].loc[:, 'MS4'])
                df_scatx[sec] = (aggsec / df_bc[name])
                df_scaty[sec] = (aggsec / aggsecemi[sec])
            elif sec == 'residential':
                df_scatx[sec] = (df_sec[name].loc[:, 'MS2']) / df_bc[name]
                aggsecemi[sec] =df_sec['emi_sec'].loc[:, 'MS2']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS2']) / (df_sec['emi_sec'].loc[:, 'MS2'])
            elif sec == 'agriculture':
                df_scatx[sec] = ((df_sec[name].loc[:, 'MS10']) / df_bc[name])
                aggsecemi[sec] = df_sec['emi_sec'].loc[:, 'MS10']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS10']) / (df_sec['emi_sec'].loc[:, 'MS10'])
            elif sec == 'transport':
                df_scatx[sec] = ((df_sec[name].loc[:, 'MS7']) / df_bc[name])
                aggsecemi[sec]=df_sec['emi_sec'].loc[:, 'MS7']
                df_scaty[sec] = (df_sec[name].loc[:, 'MS7']) / (df_sec['emi_sec'].loc[:, 'MS7'] )
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
    plt.xlabel('Potential, $\phi_p^{{Sector}}$, for {}'.format(name))
    plt.ylabel('{}, for {}'.format(yaxiskword, name))

    for ind, sec in enumerate(sect_aggr):
        patches.append(mpatches.Patch(color= 'C{}'.format(int(ind)), label='{}'.format(sec)))
    plt.legend(handles=patches, ncol=2)
    fig.savefig(path_figures + 'perprecpersectex_{}.pgf'.format(name)) #,
    fig.savefig(path_figures + 'perprecpersectex_{}.pdf'.format(name))


    ### Part 2 - areas contribution:
    yaxiskword = 'Efficiency, $\\eta_{{p,area}}$'
    marks = ['o', 'v', 's', '*','D']
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
    for indp, precursor in enumerate(prec_lst):
        marker = mlines.Line2D([], [], color='black', marker=marks[indp],
                          markersize=6, label=prec_lst[indp], linestyle = 'None')
        patches.append(marker)
        # result of regional reductions (fua, not_fua, country)
        df_scatx = pd.DataFrame(index=sources, columns=areas)
        # dataframe with y axix : efficiency (potency)
        df_scaty = pd.DataFrame(index=sources, columns=areas)
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_scaty['FUA'] = (df_res[name].ix[:, 'FUA'] /
            (df_res['emi_reg'].ix[:, 'FUA'] * precred / 100))

        df_scaty['not_FUA'] = (df_res[name].ix[:, 'not_FUA'] /
            (df_res['emi_reg'].ix[:, 'not_FUA'] * precred / 100))

        df_scatx['FUA'] = ((df_res[name].ix[:, 'FUA'] / df_bc[name]) /
                                 (precred / 100))
        df_scatx['not_FUA'] = ((df_res[name].ix[:, 'not_FUA'] / df_bc[name]) /
                                 (precred / 100))
        df_scatx.to_csv(path_figures + 'df_potential_area_{}{}.csv'.format(name, precursor))
        df_scaty.to_csv(path_figures + 'df_efficiency_area_{}{}.csv'.format(name, precursor))

        plt.plot(df_scatx['FUA'], df_scaty['FUA'], marks[indp], color= 'C{}'.format(int(1)))
        plt.plot(df_scatx['not_FUA'], df_scaty['not_FUA'], marks[indp], color= 'C{}'.format(int(2)))
    plt.xlabel('Potential, $\phi_{{p,area}}$, for {}'.format(name))
    plt.ylabel('{}, for {}'.format(yaxiskword, name))
    patches.append(mpatches.Patch(color= 'C{}'.format(int(1)), label = '{}'.format('FUA')))
    patches.append(mpatches.Patch(color= 'C{}'.format(int(2)), label = '{}'.format('non FUA')))
    plt.legend(handles=patches)

    fig.savefig(path_figures + 'perprecperareatex_{}.pgf'.format(name)) #,
    fig.savefig(path_figures + 'perprecperareatex_{}.pdf'.format(name))

def plot_ratios_fuas(name, prec_lst, sources, path_figures, path_results):
    """
    Plots the ratios of impacts for: country, fua and non fua reductions

    Created on Thu Apr 20 16:26:42 2017
    @author: peduzem
    """
    marks = ['-o', '-v', '-s', '-*', '-D']
    if not os.path.exists(path_figures):
                os.makedirs(path_figures)
    targets = ['ref_ratio', 'fua_ratio', 'not_fua_ratio']
    df_plots = pd.DataFrame(index=sources, columns=targets)
    for indp, precursor in enumerate(prec_lst):
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_plots['ref_ratio'] = (df_res[name].loc[sources].loc[:,'country'] /
            (df_res['emi_reg'].loc[sources].loc[:,'country']))

        df_plots['fua_ratio'] = (df_res[name].loc[sources].loc[:,'FUA'].where(df_res[name].loc[sources].loc[:,'FUA']>thr)/
            (df_res['emi_reg'].loc[sources].loc[:,'FUA']))

        df_plots['not_fua_ratio'] = (df_res[name].loc[sources].loc[:,'not_FUA'].where(df_res[name].loc[sources].loc[:,'not_FUA']>thr)/
            (df_res['emi_reg'].loc[sources].loc[:,'not_FUA']))
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
        p1 = plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'], marks[indp], color = '#000000')
        p2 = plt.plot(ind, df_plots['fua_ratio']/df_plots['ref_ratio'], marks[indp])
        p3 = plt.plot(ind, df_plots['not_fua_ratio']/df_plots['ref_ratio'], marks[indp])
        plt.ylabel('{} performance ratio, $\\rho_{{s,{}}}$'.format(name, precursor))
        plt.legend((p1[0], p2[0], p3[0]), ('country', 'FUA', 'Not FUA'), loc='best', ncol=1)
        plt.savefig(path_figures + 'regionalratios{}_{}.png'.format(name, precursor),
                    bbox_inches = "tight", dpi=300)
        fig.savefig(path_figures + 'regionalratiostex{}_{}.pgf'.format(name, precursor))
        fig.savefig(path_figures + 'regionalratiostex{}_{}.pdf'.format(name, precursor))
        plt.show()

#        N = len(sources)
#        ind = np.arange(N)
#        fig, ax = plt.subplots(figsize=(14, 8))
#        ax.set_ylim([0,6])
#        ax = plt.minorticks_on()
#        ax = plt.tick_params(axis='x',which='minor',bottom='off')
#        plt.xticks(ind, sources)
#        p1 = plt.plot(ind, df_plots['ref_ratio']/df_plots['ref_ratio'], '-o')
#        p2 = plt.plot(ind, df_plots['fua_ratio']/df_plots['ref_ratio'], '-o')
#        p3 = plt.plot(ind, df_plots['not_fua_ratio']/df_plots['ref_ratio'], '-o')
#        plt.ylabel('{} ratio'.format(name))
#        plt.legend((p1[0], p2[0], p3[0]), ('country', 'FUA', 'Not Fua'), loc='best', ncol=1)
#        plt.savefig(path_figures + 'regionalratios{}_{}.png'.format(name, prec),
#                    bbox_inches = "tight", dpi=300)
        plt.show()

def sectorialeffect(path_results, path_workdir, prec, precred, sources):

    if not os.path.exists(path_workdir):
        os.makedirs(path_workdir)

    targets = ['concentration', 'exposure', 'mortality', 'YOLLs',
               'emi_sec']
    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
#    ms_list = ['MS1', 'MS2', 'MS3', 'MS4', 'M7', 'M8', 'M9', 'M10']
    iterables = [sources, ms_list]
    index = pd.MultiIndex.from_product(iterables, names=['country_code', 'MS'])
    df_sec = pd.DataFrame(index=index, columns=targets)

    # Precomputation:
    # Emission per country and per sector
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()
    emi = {}
    emi_sec = {}
    for country in sources:
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        if not os.path.exists(path_areacountry_nc):
            areacountry = gridint_toarray('NUTS_Lv0', 'parea', country)
            write_nc(areacountry, path_areacountry_nc, 'AREA', '%')
        else:
            print(('WARNING: Loading previously calculated country areas'))
        (emi[country], emi_sec[country]) = precompute(path_areacountry_nc)

    # create the files path_reduction_txt with perreduction for each
    # sector and every pollutant
    for ms in ms_list:
        df_red = pd.DataFrame(np.zeros(
                              shape=(len(precursor_lst), len(ms_list))),
                              index=precursor_lst, columns=ms_list)
        path_reduction_txt = (path_workdir +
                              '\\red_ass_{}{}.txt'.format(ms, prec))
        df_red.ix[prec][ms] = precred
        df_red.to_csv(path_reduction_txt, sep='\t', index_label='POLL')

    # Computation:
    # Run model 1 with for each country and each macrosector
    # and calcluate the impacts
    for country in sources:
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        output = path_workdir + '\\output_{}_{}_{}\\'.format(prec, country,
                                                             precred)
        for it in np.arange(0, len(ms_list)):
            path_reduction_txt = path_workdir + '\\red_ass_{}{}.txt'.format(
                    ms_list[it], prec)
            outputx = output + '{}\\'.format(ms_list[it])
            if not os.path.exists(outputx):
                os.makedirs(outputx)
                proglog_filename = path_result_cdf_test + 'proglog'
                write_progress_log(proglog_filename, 25, 2)
                start = time()
                shrp.module1(path_emission_cdf_test, path_areacountry_nc,
                             path_reduction_txt, path_base_conc_cdf_test,
                             path_model_cdf_test, outputx)
        #            shrp.module1(path_emission_cdf_test, nc_redarea,
        #                         path_reduction_txt, path_model_cdf_test,
        #                         output)
                stop = time()
                print('Module 1 for {}, {}, run time: {} sec.'.format(country,
                      prec, (stop-start)))
                remove(proglog_filename)
            else:
                print('WARNING: Loading previously calculated module1 results')
            # Get results
            deltaconc = outputx + 'delta_concentration.nc'
            rootgrp = Dataset(deltaconc, mode='r')
            bc_pm25_conc = rootgrp.variables['delta_concentration'][:]
            # bc_pm25_units = rootgrp.variables['conc'].units
            rootgrp.close()
            path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
            rootgrp = Dataset(path_areacountry_nc, mode='r')
            areacountry = rootgrp.variables['AREA'][:]
            # bc_pm25_units = rootgrp.variables['conc'].units
            rootgrp.close()
            # Substitute 0s to NaN
            conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))
            # Total popluation in the country
            poptot = np.sum(popall * areacountry / 100)
            # REduction of emissions
#            emi_sec_red = emi_sec[country][prec, it] * precred / 100
            # Fill dataframe with results
            # which are moltiplied by 2 to get the total potential (as if the
            # reduction was 100%)
            df_sec.ix[country, ms_list[it]]['concentration'] = (np.sum((
                    1000 * conc * areacountry * surf / 100) /
                    np.sum(areacountry * surf / 100))) / (precred/100) #  ng/m3
            df_sec.ix[country, ms_list[it]]['exposure'] = (np.sum(( #  ng/m3
                    1000 * popall * conc * areacountry / 100) / poptot)) / (precred/100)

            (deltayll_reg, delta_mort_reg, deltayll_spec_reg) = calc_impacts(
             deltaconc, path_areacountry_nc, country, spec='d', approx='e',
             miller=True, std_life_exp=70)
            df_sec.ix[country, ms_list[it]]['mortality'] = np.asarray(
                    delta_mort_reg) / (precred/100)
            df_sec.ix[country, ms_list[it]]['YOLLs'] = np.asarray(
                    deltayll_reg) / (precred/100)
            df_sec.ix[country, ms_list[it]]['emi_sec'] = (
                    emi_sec[country][prec, it])
    df_sec.to_csv(path_results + 'df_sec_{}.csv'.format(prec))


def areaeffect(path_results, path_workdir, prec, precred, sources):
    # create work directory if it does not exist
    if not os.path.exists(path_workdir):
        os.makedirs(path_workdir)

    targets = ['concentration', 'exposure', 'mortality', 'YOLLs', 'emi_reg']
    regreductions = ['FUA', 'not_FUA', 'country']
    iterables = [sources, regreductions]
    index = pd.MultiIndex.from_product(iterables,
                                       names=['country_code', 'reg_reduction'])
    df_res = pd.DataFrame(index=index, columns=targets)

    # Precomputation:
    # create the files path_reduction_txt with perreduction for every
    # sector and every precursor
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()
    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
    for precursor in precursor_lst:
        df_red = pd.DataFrame(np.zeros(
                              shape=(len(precursor_lst), len(ms_list))),
                              index=precursor_lst, columns=ms_list)
        path_reduction_txt = path_workdir + '\\red_ass_{}.txt'.format(precursor)
        df_red.ix[precursor] = precred
        df_red.to_csv(path_reduction_txt, sep='\t', index_label='POLL')

    # Prepare areas to apply reductions
    # create multiindex dataframe adding the country code to the fuacode
    fuas = pd.read_csv('input\\selection\\fua_names.csv', index_col=[0])
    cou = [list(fuas.index.values)[fua][0:2]
           for fua in np.arange(len(list(fuas.index.values)))]
    # set country code as a column
    fuas = fuas.assign(country=cou)
    # set FUA code as another column
    fuas = fuas.assign(FUA_CODE=list(fuas.index.values))
    # set multi-index df with country and fua_code indices
    fuas = fuas.set_index(['country', 'FUA_CODE'])

    # Open model file, I just need it to load lons and lats etc.
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()

    # Create areas to apply the reduction to:
    # Area of FUAs in the country
    for country in sources:
        path_areafuas_nc = path_workdir+'\\area_fuas_{}.nc'.format(country)
        areafuas = np.zeros((1, len(lat_array), len(lon_array)))
        if not os.path.exists(path_areafuas_nc):
            for fua in fuas.ix[country].index.values:
                areafuas = areafuas + gridint_toarray('FUA_CODE', 'parea', fua)
                print(fua)
            write_nc(areafuas, path_areafuas_nc, 'AREA', '%')
            save_obj(areafuas, 'area_fuas_{}'.format(country))
        else:
            print('WARNING: Loading previously calculated fuas')
            areafuas = load_obj('area_fuas_{}'.format(country))

        # Area of the country
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        rootgrp = Dataset(path_areacountry_nc, mode='r')
        areacountry = rootgrp.variables['AREA'][:]
        # bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()
        #    areacountry = gridint_toarray('NUTS_Lv0', 'parea', country)
        #    write_nc(areacountry, path_areacountry_nc, 'AREA', '%')

        # Area of the country that does not belong to FUAs, removing numbers
        # below zero to take care of FUAS across borders (it is a problem a bit
        # for Luxemburg)
        areanonfuasraw = areacountry - areafuas
        areanonfuas = np.where((areanonfuasraw < 0), 0, (areanonfuasraw))
        path_areanonfuas_nc = (path_workdir +
                               '\\area_countrynonfuas_{}.nc'.format(country))
        write_nc(areanonfuas, path_areanonfuas_nc, 'AREA', '%')

        # Computation:
        nc_redarea = [path_areafuas_nc, path_areanonfuas_nc,
                      path_areacountry_nc]
        output = (path_workdir +
                  '\\output_{}_{}_{}\\'.format(prec, country, precred))
        path_reduction_txt = path_workdir + '\\red_ass_{}.txt'.format(prec)
        # Total emissions in each region
        emi_reg = {}
        emi_sec = {}
        for it in np.arange(0, len(nc_redarea)):
            outputx = output + '{}\\'.format(regreductions[it])
            (emi_reg[it], emi_sec[it]) = precompute(nc_redarea[it])
            # run module 1
            if not os.path.exists(outputx):
                os.makedirs(outputx)
                proglog_filename = path_result_cdf_test + 'proglog'
                write_progress_log(proglog_filename, 25, 2)
                start = time()
                shrp.module1(path_emission_cdf_test, nc_redarea[it],
                             path_reduction_txt, path_base_conc_cdf_test,
                             path_model_cdf_test, outputx)
        #       shrp.module1(path_emission_cdf_test, nc_redarea,
        #                    path_reduction_txt, path_model_cdf_test, output)
                stop = time()
                print('Module 1 for {}, {}, run time: {} sec.'.format(country,
                      prec, (stop-start)))
                remove(proglog_filename)
            else:
                print('WARNING: Loading previously calculated module1 results')

            deltaconc = outputx + 'delta_concentration.nc'
            rootgrp = Dataset(deltaconc, mode='r')
            bc_pm25_conc = rootgrp.variables['delta_concentration'][:]
            rootgrp.close()

            conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))
            poptot = np.sum(popall * areacountry / 100)

            df_res.ix[country, regreductions[it]]['concentration'] = (
                    np.sum((1000 * conc * areacountry * surf / 100) /
                           np.sum(areacountry * surf / 100))) / (precred/100)

            df_res.ix[country, regreductions[it]]['exposure'] = (
                    np.sum((1000 * popall * conc * areacountry / 100) /
                           poptot)) / (precred/100)

            (deltayll_reg, delta_mort_reg, deltayll_spec_reg) = calc_impacts(
             deltaconc, path_areacountry_nc, country, spec='d', approx='e',
             miller=True, std_life_exp=70)

            df_res.ix[country, regreductions[it]]['mortality'] = (
                    np.asarray(delta_mort_reg)) / (precred/100)
            df_res.ix[country, regreductions[it]]['YOLLs'] = (
                    np.asarray(deltayll_reg)) / (precred/100)
            df_res.ix[country, regreductions[it]]['emi_reg'] = (
                    emi_reg[it][prec])

    df_res.to_csv(path_results + 'df_res_{}.csv'.format(prec))


def basecase(path_results, path_workdir, sources):
    targets = ['concentration', 'exposure', 'mortality', 'YOLLs', 'emi']
    df_bc = pd.DataFrame(index=sources, columns=targets)

    for country in sources:
        # Area of the country
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        rootgrp = Dataset(path_areacountry_nc, mode='r')
        areacountry = rootgrp.variables['AREA'][:]
        path_bc_conc_nc = path_base_conc_cdf_test
        rootgrp = Dataset(path_bc_conc_nc, mode='r')
        bc_pm25_conc = rootgrp.variables['conc'][:]
        # bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()
        conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))
        poptot = np.sum(popall * areacountry / 100)
        df_bc.ix[country]['concentration'] = (np.sum((
                1000 * conc * areacountry * surf / 100) /
                np.sum(areacountry * surf / 100)))
        df_bc.ix[country]['exposure'] = (np.sum((
                1000 * popall * conc * areacountry / 100) / poptot))
        (deltayll_reg, delta_mort_reg, deltayll_spec_reg) = calc_impacts(
         path_bc_conc_nc, path_areacountry_nc, country, spec='a', approx='e',
         miller=True, std_life_exp=70)
        df_bc.ix[country]['mortality'] = np.asarray(
                delta_mort_reg)
        df_bc.ix[country]['YOLLs'] = np.asarray(
                deltayll_reg)

    df_bc.to_csv(path_results + 'df_bc.csv')
# -------------------------------------------------------------------------
# main program starts here
# -------------------------------------------------------------------------

if __name__ == '__main__':

    # -------------------------------------------------------------------------
    # Input data
    # -------------------------------------------------------------------------

    # Population to calculate exposure:
    path_tiff = 'input/pop/7km_Qi_2010.tif'
    popall = tiftogridgeneral(path_tiff)
    path_surf_nc = 'input/JRC01.nc'

    # Area of each cell in the domain:
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    # area_units = rootgrp.variables['surface'].units
    rootgrp.close()

    # -------------------------------------------------------------------------
    # General inputs
    # -------------------------------------------------------------------------
#    sources = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES',
#               'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
#               'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
    sources = ['BE', 'DE','ES', 'FR', 'IT', 'UK']
    prec_lst = ['NOx', 'NMVOC', 'NH3', 'PPM', 'SOx']
#    prec = 'PPM'  # Precursor
    precred = 50  # Percentage reduction of precursors

    # -------------------------------------------------------------------------
    # 1 - SECTORIAL EFFECT
    # -------------------------------------------------------------------------
    # Prepare structure to save results
    path_results = 'spaflextest\\'
    path_workdir = 'spaflextest\\workdirsec'
    for prec in prec_lst:
        sectorialeffect(path_results, path_workdir, prec, precred, sources)
    # -------------------------------------------------------------------------
    # 2 - REGIONAL EFFECT
    # -------------------------------------------------------------------------
    # Prepare structure to save results
    path_results = 'spaflextest\\'
    path_workdir = 'spaflextest\\workdirreg'
    for prec in prec_lst:
        areaeffect(path_results, path_workdir, prec, precred, sources)

    # -------------------------------------------------------------------------
    # 4 - Base case for comparison
    # -------------------------------------------------------------------------
    basecase(path_results, path_workdir, sources)


    # -------------------------------------------------------------------------
    # 5 - Display results
    # -------------------------------------------------------------------------

#    sources = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES',
#           'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
#           'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']
#    sources = ['AT', 'BE', 'DE', 'ES', 'FR', 'IT', 'UK']
    thr= 50 #threshold - minimum value of absolute potential in ng/m3 to include resutls
    path_figures = path_results
#    path_results = 'spaflextest\\'
#    prec = 'PPM'
    N = len(sources)
    ind = np.arange(N)

    targets = ['ref_ratio']
    # names = ['concentration', 'exposure']
    names = ['exposure', 'concentration']
    # ms_list = ['MS1', 'MS2', 'MS7', 'MS10']
    # PLOT RATIOS for MacroSectors

    df_bc = pd.read_csv(path_results + 'df_bc.csv', index_col=[0])
    for name in names:
        plot_ratios_fuas(name, prec_lst, sources, path_figures, path_results)
    for name in names:
        plot_ratios_secs(name, prec_lst, sources, path_figures, path_results)
    for name in names:
        plot_ratios_aggsec(name, prec_lst, sources, path_figures, path_results)
    for name in names:
        plot_effpot(prec_lst, name, sources, path_results,
                  path_figures)

    # Reviewed until here


    # PLOT FUAs and Non-FUAs contributions
    for prec in prec_lst:
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(prec), index_col=[0, 1])
        for name in names:
            datafua = df_res[name].ix[:, 'FUA']
            datanfua = df_res[name].ix[:, 'not_FUA']
            ind = np.arange(N)    # the x locations for the groups
            width = 0.35    # the width of the bars: can also be len(x) sequence
#            fig, ax = plt.subplots(figsize=(14, 8))
            fig= plt.figure(figsize=figsize(1))
            ax = fig.add_subplot(111)
            ax = plt.minorticks_on()
            ax = plt.tick_params(axis='x',which='minor',bottom='off')
            p1 = plt.bar(ind, datafua, width, color='#d62728')
            p2 = plt.bar(ind, datanfua, width,
                         bottom = datafua)
            plt.ylabel('{} [ng/m3]'.format(name))
            plt.title('FUAs and Not-FUAs contribution')
            plt.xticks(ind, sources)
            # lt.yticks(np.arange(0, 81, 10))
            plt.legend((p1[0], p2[0]), ('FUA', 'Not Fua'), loc='best', ncol=1)  #bbox_to_anchor=(1.13, 0.6),
            plt.savefig(path_results + 'regional{}_{}.png'.format(name, prec),
                        bbox_inches = "tight", dpi=300)
            plt.show()

    # plot MS contributions
    for prec in prec_lst:
        df_sec = pd.read_csv(path_results +
                        'df_sec_{}.csv'.format(prec), index_col=[0, 1])
        for name in names:
            width = 0.35    # the width of the bars: can also be len(x) sequence
#            fig, ax = plt.subplots(figsize=(14, 8))
            fig= plt.figure(figsize=figsize(1))
            ax = fig.add_subplot(111)
            ax = plt.minorticks_on()
            ax = plt.tick_params(axis='x',which='minor',bottom='off')
            epms = df_sec[name]
            bot = 0
            for ms in np.arange(1, 11):
                msdata = epms.ix[:, 'MS{}'.format(ms)]
                p = plt.bar(ind, (msdata ),
                            width, bottom=bot, label='MS{}'.format(ms))
        #        plt.legend((p[0]), ('MS{}'.format(ms)))!
                bot = bot + (msdata)
            plt.legend(ncol=2, loc='best') # , bbox_to_anchor=(1.13, 0.6)
            plt.ylabel('{} [ng/m3]'.format(name))
            plt.title('Sectors contribution')
            plt.xticks(ind, sources)
            plt.savefig(path_results + 'sectors{}_{}.png'.format(name, prec),
                        bbox_inches = "tight", dpi=300)
            plt.show()

    # Plot contribution of precursors:
    for prec in prec_lst:
        for name in names:  # TODO check
            width = 0.35    # the width of the bars: can also be len(x) sequence
#            fig, ax = plt.subplots(figsize=(14, 8))
            fig= plt.figure(figsize=figsize(1))
            ax = fig.add_subplot(111)
            ax = plt.minorticks_on()
            ax = plt.tick_params(axis='x',which='minor',bottom='off')
            epms = df_bc[name]
            bot = df_res[name].ix[:, 'country']
            p2 = plt.bar(ind, epms, width, bottom = 0)
            p1 = plt.bar(ind, bot, width)
            plt.ylabel('{} [ng/m3]'.format(name))
            plt.title('Precursor contribution')
            plt.xticks(ind, sources)
            plt.legend((p1[0], p2[0]), ('{} contribution'.format(prec), 'Other contribution'), loc='best', ncol=1) # bbox_to_anchor=(1.01, 0.6),
            plt.savefig(path_results + 'basecase{}_{}.png'.format(name,prec),
                        bbox_inches = "tight", dpi=300)
            plt.show()

    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()
    # for name in names:
    name = 'concentration'
#    ms_list = ['MS1', 'MS2', 'MS7', 'MS10']#, 'MS2', 'MS7', 'MS10']
#    name = 'concentration'
    precursor_lst = ['PPM', 'SOx', 'NOx', 'NH3']#, 'NMVOC']
    marks = ['o', 'v', 's', '*']
    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
    # Plots by sectors:
#    fig, ax = plt.subplots(figsize=(14, 8))
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax = plt.minorticks_on()
    plt.ylim([0,0.205])
    plt.xlim([0,8])
    patches = []
    for indp, precursor in enumerate(precursor_lst):
        df_bc = pd.read_csv(path_results +
                        'df_bc.csv', index_col=[0])
        df_res = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_sec = pd.read_csv(path_results +
                        'df_sec_{}.csv'.format(precursor), index_col=[0, 1])
        df_plots =  pd.read_csv(path_results +
                        'df_plots_{}{}.csv'.format(name, precursor),
                                      index_col=[0])
        df_scatx = pd.DataFrame(index=sources, columns=ms_list)
        df_scaty = pd.DataFrame(index=sources, columns=ms_list)
        marker = mlines.Line2D([], [], color='black', marker=marks[indp],
                          markersize=6, label=precursor_lst[indp], linestyle = 'None')
        patches.append(marker)
        for ind, ms in enumerate(ms_list):
            # df_scaty[ms] = df_sec[name].ix[:, ms]/df_res[name].ix[:, 'country']
            df_scaty[ms] = (df_sec[name].ix[:, ms])/df_bc[name]
            df_scatx[ms] =  df_plots['{}'.format(ms)]/df_plots['ref_ratio']
            plt.plot(df_scatx[ms], df_scaty[ms], marks[indp], label = '{}'.format(ms), color= 'C{}'.format(int(ms[2:])-1))
#    plt.ylabel('"potential" $\Phi_p^{{MS}}$ = $\Delta$ {}$_{{p,MS}}$/{}{}'.format(effect, effect, '$^{bc}$'))
    plt.ylabel('Potential, $\phi_p^{{MS}}$, for {}'.format(name))
    plt.xlabel('Efficiency ratio, $\\rho_p^{{MS}}$, for {}'.format(name))
    for ind, ms in enumerate(ms_list):
            patches.append(mpatches.Patch(color= 'C{}'.format(int(ms[2:])-1), label='{}'.format(ms)))
    plt.legend(handles=patches)
    plt.savefig(path_results + 'perprecpersec_{}.png'.format(name),
                    bbox_inches = "tight", dpi=300)
    # plt.legend(markerscale=0.7, scatterpoints=1)
    plt.show()

    df_tot_exp_ms = pd.DataFrame(np.zeros((len(sources),len(ms_list))),
                                 index=sources, columns=ms_list)
    df_tot_exp_prec = pd.DataFrame(np.zeros((len(sources),len(precursor_lst))),
                                 index=sources, columns=precursor_lst)
    exp_prec_ratio = pd.DataFrame(np.zeros((len(sources),len(precursor_lst))),
                                 index=sources, columns=precursor_lst)
    exp_ms_ratio = pd.DataFrame(np.zeros((len(sources),len(ms_list))),
                                 index=sources, columns=ms_list)
    df_bc = pd.read_csv(path_results + 'df_bc.csv', index_col=[0])

    for ind, precursor in enumerate(precursor_lst):

        df_res = pd.read_csv(path_results +
                         'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_sec = pd.read_csv(path_results +
                         'df_sec_{}.csv'.format(precursor), index_col=[0, 1])
        df_plots =  pd.read_csv(path_results +
                            'df_plots_{}{}.csv'.format(name, precursor),
                                      index_col=[0])
#        df_tot_exp_prec[precursor] = df_sec[name].ix[:, ms]
        df_tot_exp_prec[precursor] = df_res[name].ix[:, 'country'] + df_tot_exp_prec[precursor]
        for ms in ms_list:
#        ms = 'MS1'
                # df_scaty[ms] = df_sec[name].ix[:, ms]/df_res[name].ix[:, 'country']
            df_tot_exp_ms[ms] = df_sec[name].ix[:, ms]+df_tot_exp_ms[ms]

    df_tot_exp_prec['sum'] = df_tot_exp_prec.sum(axis=1)
    df_tot_exp_ms['sum'] = df_tot_exp_ms.sum(axis=1)
    ### plots
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax = plt.minorticks_on()
    for ind, precursor in enumerate(precursor_lst):
        exp_prec_ratio[precursor] = df_tot_exp_prec[precursor] / df_tot_exp_prec['sum']
        plt.plot((np.ones(len(sources))*ind), exp_prec_ratio[precursor],
                 marks[ind], label = '{}'.format(precursor), color = 'black')
    plt.xticks(np.arange(len(precursor_lst)), precursor_lst)
    plt.ylabel('$\sum_{{MS}}${}/$\sum_{{MS,p}}${}'.format('Ex', 'Ex'))
    plt.xlabel('Precursor')
    plt.savefig(path_results + 'preccontr_{}.png'.format(name),
                    bbox_inches = "tight", dpi=300)
    plt.legend()
    plt.show()

    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax = plt.minorticks_on()
    for ind, ms in enumerate(ms_list):
        exp_ms_ratio[ms]=df_tot_exp_ms[ms] / df_tot_exp_ms['sum']
        fig2 = plt.plot((np.ones(len(sources))*ind), exp_ms_ratio[ms], 'o', label = '{}'.format(ms))
    plt.xticks(np.arange(len(ms_list)), ms_list)
    plt.ylabel('$\sum_{{p}}${}/$\sum_{{MS,p}}${}'.format('Ex', 'Ex'))
    plt.xlabel('Macro Sectors')
    plt.savefig(path_results + 'macroseccontr_{}.png'.format(name),
                    bbox_inches = "tight", dpi=300)
    plt.legend(ncol=2, loc='best')
    plt.show()


