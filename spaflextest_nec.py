# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 15:28:00 2017
Scipt to test the national emissions reduction commitments
Last version! 

@author: peduzem
"""

import matplotlib as mpl
mpl.use('pgf')

# import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap

from netCDF4 import Dataset
import numpy as np
import os as os
from os import remove
import pandas as pd  # conda install pandas
from pandas import ExcelWriter
from time import time  # for module1

from sherpa_auxiliaries import (write_progress_log)
from sherpa_globals import (path_base_conc_cdf_test, path_emission_cdf_test,
                            path_model_cdf_test, path_result_cdf_test)
from sherpa_auxiliaries_epe import (gridint_toarray, write_nc,
                                    tiftogridgeneral)
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
import module1 as shrp

def read_inventory(year, em_inv, precursors):
    """
    Created on Wed Jun  7 15:28:00 2017
    Script to read the inventories of TNO EMEP or GAINS
    - input: year, em_inv, precursors
    - supporting files: excel files with the emission inventories
    - output: dataframe with emissions per country as reported in the inventory
    @author: peduzem
    """

    # nomenclature in inventory : EMEP and TNO
    sourcesinv = ['AT', 'BE', 'BG', 'CY', 'CZ', 'DE', 'DK', 'EE', 'GR', 'ES',
                  'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'MT',
                  'NL', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'GB']

   # Create a dataframe with the emissions per precursors for the reference
   # year

    # Dataframe with emissions per country as reported in the inventory
    df_e_inv_eu = pd.DataFrame(index=sourcesinv, columns=precursors)

    # Build df_e_inv_eu by reading the inventory file
    if em_inv == 'GAINS':
        path_emi_inv = 'input\\gains_inventory.xlsx'
        rename_dct = pd.read_excel('input\\utilities\\country_codes_4to2.xlsx',
                                   index_col=0, header=None).to_dict()
        for precu in precursors:
            if precu == 'PPM':
                s_name = 'p25'
            if precu == 'NOx':
                s_name = 'nox'
            if precu == 'SOx':
                s_name = 'so2'
            if precu == 'NH3':
                s_name = 'nh3'
            if precu == 'NMVOC':
                s_name = 'voc'
            df_e_inv = pd.read_excel(path_emi_inv,
                                     sheetname='{}'.format(s_name),
                                     skiprows=6,
                                     header=0, skip_footer=44,
                                     index_col=[0, 2])
            df_e_inv.rename(index=rename_dct[1], inplace=True)
            df_e_inv_eu[precu] = df_e_inv[year].loc[:, 'SUM']
    elif em_inv == 'EMEP' or em_inv == 'TNO':
        if em_inv == 'EMEP':
            path_emi_inv = 'input\\emep_emissions.xlsx'
        elif em_inv == 'TNO':
            path_emi_inv = 'input\\tno_emissions.xlsx'
        for precu in precursors:
            if precu == 'PPM':
                df_e_inv = pd.read_excel(path_emi_inv, sheetname='PM2.5',
                                         index_col=[1])
                df_e_inv_eu[precu] = df_e_inv.ix[sourcesinv][year]
            else:
                df_e_inv = pd.read_excel(path_emi_inv,
                                         sheetname='{}'.format(precu),
                                         index_col=[1])
                df_e_inv_eu[precu] = df_e_inv.ix[sourcesinv][year]

    # inventory for EU28:
    df_e_inv_eu.ix['EU28'] = df_e_inv_eu.ix[sourcesinv].sum(axis=0)
    
    # Rename keys so they are the same as in SHERPA..
    df_e_inv_eu.rename(index={'GR': 'EL', 'GB': 'UK'}, inplace=True)
    
    return df_e_inv_eu

def nec_app(an_ty, opt, name, sources, path_results, path_figures, em_inv):
    '''Ideal application of the NEC directive.
       INPUT:
           - an_ty: analysis type sec or spa
           - opt: I or II
           - name: exposure 
           - sources: country codes of interests 
           - path_results:
           - path_figures
           - em_inv: emission inventory name (e.g. GAINS)
       OUTPUT
    '''

    # Prepare dictionaries for country-to-grid data
    df_res = {}

    # dataframe with the efficiencies for each macrosector/area and precursor
    df_eff = {}
    df_pot = {}
    df_eff_country = {}
    df_pot_country = {}

    # Relative efficiency with respect to the maximum one
    df_releff = {}
    # relative potential with respect to the maximum one
    df_relpot = {}
    # Dataframe with the indicators
    df_ind = {}
    df_emi = {}

    sect_aggr = ['industry', 'residential', 'agriculture', 'transport',
                 'other']
    areas = ['FUA', 'not_FUA']

    if an_ty == 'sec':
        columnsresults = sect_aggr
    elif an_ty == 'spa':
        columnsresults = areas
        
    # Grid to grid data:
    # reading data and preparing dataframe for results:
    for precursor in precursors:
        if an_ty == 'sec':
            df_eff[precursor] = pd.read_csv(
                   path_figures + 'df_efficiency_aggr_{}{}.csv'.format(
                    name, precursor), index_col=[0])
            df_pot[precursor] = pd.read_csv(path_figures + 'df_potential_aggr_{}{}.csv'.format(name, precursor), index_col=[0])
            df_emi[precursor] = pd.read_csv(path_figures + 'df_aggsecemi_aggr_{}.csv'.format(precursor), index_col=[0])
        elif an_ty == 'spa':
            df_eff[precursor] = pd.read_csv(path_figures + 'df_efficiency_area_{}{}.csv'.format(name, precursor), index_col=[0])
            df_pot[precursor] = pd.read_csv(path_figures + 'df_potential_area_{}{}.csv'.format(name, precursor), index_col=[0])
            df_emi[precursor] = pd.read_csv(path_results + 'df_res_{}.csv'.format(precursor), index_col=[0, 1], usecols = [0,1,6])

        df_releff[precursor] = pd.DataFrame(index=list(df_eff[precursor].index),
                                 columns = columnsresults)
        df_relpot[precursor] = pd.DataFrame(index=list(df_pot[precursor].index),
                                 columns = columnsresults)
        df_ind[precursor] = pd.DataFrame(index=list(df_pot[precursor].index),
                                 columns = columnsresults)
   
    # Calculating the performance indicator
        for sect in columnsresults:
            df_releff[precursor][sect] = (df_eff[precursor]['{}'.format(sect)] /
                                             df_eff[precursor][columnsresults].max(axis=1))
            df_relpot[precursor][sect] = (df_pot[precursor]['{}'.format(sect)] /
                                             df_pot[precursor][columnsresults].max(axis=1))
            if opt == 'I':
                df_ind[precursor][sect] = df_releff[precursor][sect] #0.5*df_relpot[precursor][sect] + 0.5 *
            elif opt == 'II':
                df_ind[precursor][sect] = df_relpot[precursor][sect]
        if opt == 'I':
            df_releff[precursor]['max'] = df_eff[precursor][columnsresults].max(axis=1)*1000
        if opt == 'II':
            df_relpot[precursor]['max'] = df_pot[precursor][columnsresults].max(axis=1)*100

    df_bc = pd.read_csv(path_results + 'df_bc.csv', index_col=[0])

    # Country to grid data:
    for precursor in precursors:
        df_res[precursor] = pd.read_csv(path_results +
                        'df_res_{}.csv'.format(precursor), index_col=[0, 1])
        df_eff_country[precursor]=(df_res[precursor][name].ix[sources].ix[:,'country'] /
            (df_res[precursor]['emi_reg'].ix[sources].ix[:,'country']))
        df_pot_country[precursor] = pd.DataFrame(index=list(df_pot[precursor].index),
                                 columns = columnsresults)

        for sect in columnsresults:
            if an_ty == 'sec':
                df_pot_country[precursor][sect] = df_eff_country[precursor] * df_emi[precursor][sect] /  df_bc[name].ix[sources]
            if an_ty == 'spa':
                df_pot_country[precursor][sect] = df_eff_country[precursor] * df_emi[precursor]['emi_reg'].loc[:,sect].loc[sources] / df_bc[name].ix[sources]
    df_tot_imp_prec = pd.DataFrame(np.zeros((len(df_pot[precursor].index),len(precursors))),
                                 index=df_pot[precursor].index, columns=precursors)
#    exp_prec_ratio = pd.DataFrame(np.zeros((len(df_pot[precursor].index),len(precursors))),
#                                 index=df_pot[precursor].index, columns=precursors)
    for ind, precursor in enumerate(precursors):
        df_pot[precursor]['sum'] = df_pot[precursor][columnsresults].sum(axis=1)
        df_tot_imp_prec[precursor] = df_pot[precursor][columnsresults].sum(axis=1)

    # prioritizing the sectors
    df_imp_de = {}
    df_imp_as = {}
    df_imp_co = {}
    lst ={}
    lst2 = {}
    red = {}
    redpertable = {}
    results = {}
    redms = {} # % reduction per macrosector
    redms2 = {}
    colres = ['case_de', 'case_as', 'countryred']

    for precursor in precursors:
        df_imp_de[precursor] = pd.DataFrame(index=[sources], columns = columnsresults)
        df_imp_as[precursor] = pd.DataFrame(index=[sources], columns = columnsresults)
        df_imp_co[precursor] = pd.DataFrame(index=[sources], columns = columnsresults)
        results[precursor] = pd.DataFrame(index=[sources], columns = colres)
        redms[precursor] = pd.DataFrame(index=[sources], columns = columnsresults)
        redms2[precursor] = pd.DataFrame(index=[sources], columns = columnsresults)
#    em_inv = 'GAINS' # 'EMEP', 'TNO', 'GAINS'
    df_e_inv_eu = read_inventory(year, em_inv, precursors)
    for precursor in precursors:
        # Calculate the absolute value of emissions after implementation of NEC:
        df_e_red_eu[precursor] = (1 - df_nec[yearred][precursor].ix[sources]) * df_e_inv_eu[precursor].ix[sources]*1000
        # Absolute value of emission reduction from year 2010
        if an_ty == 'sec':
            if em_inv == 'GAINS':
                emi2010 = (read_inventory(2010, em_inv, [precursor])*1000).loc[sources][precursor]
            else:
                emi2010 = df_emi[precursor].sum(axis=1)
            red[precursor] = emi2010 - df_e_red_eu[precursor]
#            redper[precursor] =red[precursor]/df_emi[precursor].sum(axis=1) * 100
            redpertable[precursor] = red[precursor]/emi2010 * 100 
        elif an_ty == 'spa':
            if em_inv == 'GAINS':
                emi2010 = (read_inventory(2010, em_inv, [precursor])*1000).loc[sources][precursor]
            else:
                emi2010 = df_emi[precursor]['emi_reg'].loc[:, 'country'].loc[sources]
            red[precursor] = emi2010.loc[sources] - df_e_red_eu[precursor]
#            redper[precursor] = red[precursor] / df_emi[precursor]['emi_reg'].loc[:, 'country'].loc[sources] * 100
            redpertable[precursor] = red[precursor]/emi2010 * 100

    for precursor in precursors:
        for country in sources:
            # order of sectors to which we can apply the reduction
            lst[country] = df_ind[precursor].ix[country].sort_values(ascending=False)
            rawvalue = red[precursor].loc[country]
            for ordms in lst[country].index:
                if an_ty == 'sec':
                    emi = df_emi[precursor][ordms].loc[country]
                elif an_ty == 'spa':
                    emi = df_emi[precursor]['emi_reg'].loc[:, ordms].loc[country]
                if emi == 0:
                    redms[precursor][ordms].loc[country] = 0
                else:
                    if rawvalue > emi:
                       redms[precursor][ordms].loc[country] = 100
                       rawvalue = rawvalue -emi
                       pass
                    elif rawvalue < emi:
                       redms[precursor][ordms].loc[country] = rawvalue / emi * 100
                       rawvalue = 0
                       pass
                    elif rawvalue == 0:
                       redms[precursor][ordms].loc[country] = 0
                       rawvalue = 0
                       pass
    for precursor in precursors:
        for country in sources:
            lst2[country] = df_ind[precursor].ix[country].sort_values(ascending=True)
            rawvalue = red[precursor].loc[country]
            for ordms in lst2[country].index:
                if an_ty == 'sec':
                    emi = df_emi[precursor][ordms].loc[country]
                elif an_ty == 'spa':
                    emi = df_emi[precursor]['emi_reg'].loc[:, ordms].loc[country]

                if emi == 0:
                    redms2[precursor][ordms].loc[country] = 0
                else:
                    if rawvalue > emi:
                       redms2[precursor][ordms].loc[country] = 100
                       rawvalue = rawvalue - emi
                       pass
                    elif rawvalue < emi:
                       redms2[precursor][ordms].loc[country] = rawvalue / emi * 100
                       rawvalue = 0
                       pass
                    elif rawvalue == 0:
                       redms2[precursor][ordms].loc[country] = 0
                       rawvalue = 0
                       pass
    for precursor in precursors:
        for sect in columnsresults:
            if an_ty == 'sec':
                emi = df_emi[precursor].ix[:, sect]
            if an_ty == 'spa':
                emi = df_emi[precursor]['emi_reg'].loc[:, sect].loc[sources]

            df_imp_de[precursor][sect] = (df_eff[precursor].ix[:, sect] * ((100 - redms[precursor][sect])/100) * emi)
            df_imp_as[precursor][sect] = (df_eff[precursor].ix[:, sect] * ((100 - redms2[precursor][sect])/100) * emi)

    for precursor in precursors:
    # total impact in the three cases per precursor
        df_imp_de[precursor]['sum']=df_imp_de[precursor][columnsresults].sum(axis=1)
        df_imp_as[precursor]['sum']=df_imp_as[precursor][columnsresults].sum(axis=1)
        df_imp_co[precursor] = df_bc[name].loc[sources]*df_tot_imp_prec[precursor]- df_eff_country[precursor] * red[precursor]

    for precursor in precursors:
    # results show the percentage reduction with respect to the base case that is obtained in each country for each precursor
        results[precursor]['case_de'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_de[precursor]['sum'])/df_bc[name])*100
        results[precursor]['case_as'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_as[precursor]['sum'])/df_bc[name])*100
        results[precursor]['countryred'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_co[precursor])/df_bc[name])*100

    results['sum'] = pd.DataFrame(np.zeros((len(sources),len(colres))), index=[sources], columns = colres)
    for col in results[precursors[0]].columns:
        for precursor in precursors:
            for country in sources:
                results['sum'][col].loc[country] = results['sum'][col].loc[country] + results[precursor][col].loc[country]
#    results[precursors[0]].columns

    for key in results.keys():
        results[key]['ratio_de'] = results[key]['case_de'] / results[key]['countryred']
        results[key]['ratio_as'] = results[key]['case_as'] / results[key]['countryred']
        results[key]['ratio_co'] = results[key]['countryred'] / results[key]['countryred']
#        results[precursor]['sectorial_de'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_de[precursor]['sum'])/(df_bc[name].loc[sources]*df_tot_imp_prec[precursor]))*100
#        results[precursor]['sectorial_as'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_as[precursor]['sum'])/(df_bc[name].loc[sources]*df_tot_imp_prec[precursor]))*100
#        results[precursor]['countryred'] = ((df_bc[name]*df_tot_imp_prec[precursor]- df_imp_co[precursor])/(df_bc[name].loc[sources]*df_tot_imp_prec[precursor]))*100

    # Display of results
#    marks = ['-o', '-v', '-s', '-*']
    N = len(sources)
    ind = np.arange(N)
    if opt == 'I':
        indicator = '\\eta'
    elif opt == 'II':
        indicator = '\\phi'

    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()
#    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_ylim([0,2])
#    ax = plt.minorticks_on()
    ax.tick_params(axis='x',which='minor',bottom='off')
    plt.xticks(ind, sources)
    plt.plot(ind, results['sum']['ratio_co'], "-" , color = 'black', label = 'flat red')
    plt.plot(ind, results['sum']['ratio_de'], "D" , color = 'purple', label = u'decreasing ${}_{{p,s}}$'.format(indicator))
    plt.plot(ind, results['sum']['ratio_as'], "D" , color = 'dodgerblue', label = u'increasing ${}_{{p,s}}$'.format(indicator))
    plt.legend(ncol=1, loc='best')
    plt.ylabel('{} performance ratio, $\\rho_{{NEC}}$'.format(name))
    fig.savefig(path_figures + 'NEC{}{}_ratio{}.png'.format(name, an_ty, opt),
                bbox_inches = "tight", dpi=300)
    fig.savefig(path_figures + 'NEC{}{}_ratio{}.pgf'.format(name, an_ty, opt))
    fig.savefig(path_figures + 'NEC{}{}_ratio{}.pdf'.format(name, an_ty, opt))
    plt.show()

    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()
#        fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_ylim([0,50])
#        ax = plt.minorticks_on()
    ax.tick_params(axis='x',which='minor',bottom='off')
    plt.xticks(ind, sources)
    plt.plot(ind, results['sum']['countryred'], "-" , color = 'black', label = 'flat red')
    plt.plot(ind, results['sum']['case_de'], "D" , color = 'purple', label = u'decreasing ${}_{{p,s}}$'.format(indicator))
    plt.plot(ind, results['sum']['case_as'], "D" , color = 'dodgerblue', label = u'increasing ${}_{{p,s}}$'.format(indicator))
    plt.legend(ncol=1, loc='best')
    plt.ylabel(u'{} performance, $\\Delta I_{{NEC}}/I$ \%'.format(name))
    fig.savefig(path_figures + 'NEC{}{}{}.png'.format(name, an_ty, opt),
                bbox_inches = "tight", dpi=300)
    fig.savefig(path_figures + 'NEC{}{}{}.pgf'.format(name, an_ty, opt))
    fig.savefig(path_figures + 'NEC{}{}{}.pdf'.format(name, an_ty, opt))
    plt.show()
    # Printing tables:
    save_xls(redms, (path_figures + 'redms_de{}{}{}.xls'.format(name, an_ty, opt)))
    save_xls(redms2, (path_figures + 'redms_as{}{}{}.xls'.format(name, an_ty, opt)))

    
    
    
def perexppercountry(name, an_ty, opt, country, order):
    """ calculate mapped percentage exposure reduction for each option"""
    # Write reduction text
    sect_aggr = ['industry', 'residential', 'agriculture', 'transport', 'other']
    areatypes = ['FUA','not_FUA']
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()

    ms_list = ['MS{}'.format(snap) for snap in np.arange(1, 11)]
    df_red = pd.DataFrame(np.zeros(
                              shape=(len(precursor_lst), len(ms_list))),
                              index=precursor_lst, columns=ms_list)
    if an_ty == 'sec':
        df_red_agg = pd.DataFrame(index=precursors, columns=sect_aggr)
    if an_ty == 'spa':
         df_red_agg = pd.DataFrame(index=precursors, columns=areatypes)
    for prec in precursors:
        if order =='dec':
            df_prec = pd.read_excel(path_figures + 'redms_de{}{}{}.xls'.format(name, an_ty, opt), prec, index_col=0)
            df_red_agg.loc[prec]=df_prec.loc[country]
        if order == 'inc':
            df_prec = pd.read_excel(path_figures + 'redms_as{}{}{}.xls'.format(name, an_ty, opt), prec, index_col=0)
            df_red_agg.loc[prec]=df_prec.loc[country]



    if an_ty == 'sec':
        path_reduction_txt = path_results + '\\redms_{}{}{}{}{}.txt'.format(country,name,order,an_ty,opt)
        for prec in precursor_lst:
#            if prec == 'NMVOC':
#                df_red.loc[prec] = 0 * np.arange(len(ms_list))
#            else:
            df_red['MS1'].loc[prec] = df_red_agg['industry'].loc[prec]
            df_red['MS3'].loc[prec] = df_red_agg['industry'].loc[prec]
            df_red['MS4'].loc[prec] = df_red_agg['industry'].loc[prec]
            df_red['MS2'].loc[prec] = df_red_agg['residential'].loc[prec]
            df_red['MS7'].loc[prec] = df_red_agg['transport'].loc[prec]
            df_red['MS10'].loc[prec] = df_red_agg['agriculture'].loc[prec]
            df_red['MS5'].loc[prec] = df_red_agg['other'].loc[prec]
            df_red['MS6'].loc[prec] = df_red_agg['other'].loc[prec]
            df_red['MS8'].loc[prec] = df_red_agg['other'].loc[prec]
            df_red['MS9'].loc[prec] = df_red_agg['other'].loc[prec]


        df_red.to_csv(path_reduction_txt, sep='\t', index_label='POLL')
        
        # Computation:
        # Run model 1 with for each country and each macrosector
        # and calcluate the impacts
        # for country in sources:
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        outputfin = path_results + 'NEC_{}{}{}{}{}\\'.format(country,name,order,an_ty,opt)
        if not os.path.exists(outputfin):
            os.makedirs(outputfin)
            proglog_filename = path_result_cdf_test + 'proglog'
            write_progress_log(proglog_filename, 25, 2)
            start = time()
            shrp.module1(path_emission_cdf_test, path_areacountry_nc,
                                     path_reduction_txt, path_base_conc_cdf_test,
                                     path_model_cdf_test, outputfin)
                #            shrp.module1(path_emission_cdf_test, nc_redarea,
                #                         path_reduction_txt, path_model_cdf_test,
                #                         output)
            stop = time()
            print('Module 1 for {}, {}, run time: {} sec.'.format(country, prec,
                  (stop-start)))
            remove(proglog_filename)
        else:
            print('WARNING: previously calculated results will be used')

        # Get results
        deltaconc = outputfin + 'delta_concentration.nc'
        rootgrp = Dataset(deltaconc, mode='r')
        d_pm25_conc = rootgrp.variables['delta_concentration'][:]
        rootgrp.close()
        path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
        rootgrp = Dataset(path_areacountry_nc, mode='r')
    #    areacountry = rootgrp.variables['AREA'][:]
        # bc_pm25_units = rootgrp.variables['conc'].units
        rootgrp.close()
        # Substitute 0s to NaN
        conc = np.where(np.isnan(d_pm25_conc), 0, (d_pm25_conc))

    if an_ty == 'spa':
        print(country)
        conc = np.zeros((1, len(lat_array), len(lon_array)))
        outputfin = path_results + 'NEC_{}{}{}{}{}\\'.format(country,name,order,an_ty,opt)
        if not os.path.exists(outputfin):
            os.makedirs(outputfin)
        for areatype in areatypes:
            output = path_results + 'NEC_{}{}{}{}{}{}\\'.format(country,name,order,an_ty,opt,areatype)
            path_reduction_txt = path_results + '\\redms_{}{}{}{}{}{}.txt'.format(country,name,order,an_ty,opt,areatype)
            for prec in precursor_lst:
#                if prec == 'NMVOC':
#                    df_red.loc[prec] = 0 * np.arange(len(ms_list))
#                else:
                df_red['MS1'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS3'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS4'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS2'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS7'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS10'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS5'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS6'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS8'].loc[prec] = df_red_agg[areatype].loc[prec]
                df_red['MS9'].loc[prec] = df_red_agg[areatype].loc[prec]

            df_red.to_csv(path_reduction_txt, sep='\t', index_label='POLL')
            # Computation:
            # Run model 1 with for each country and each macrosector
            # and calcluate the impacts
            # for country in sources:
            if areatype == 'FUA':
                path_fua_country = 'D:\\sherpa.git\\Sherpa\\asstest\\workdirreg\\'
                path_areacountry_nc = path_fua_country + 'area_fuas_{}.nc'.format(country)
            elif areatype == 'not_FUA':
                path_fua_country = 'D:\\sherpa.git\\Sherpa\\asstest\\workdirreg\\'
                path_areacountry_nc = path_fua_country + 'area_countrynonfuas_{}.nc'.format(country)

            if not os.path.exists(output):
                os.makedirs(output)
                proglog_filename = path_result_cdf_test + 'proglog'
                write_progress_log(proglog_filename, 25, 2)
                start = time()
                shrp.module1(path_emission_cdf_test, path_areacountry_nc,
                                         path_reduction_txt, path_base_conc_cdf_test,
                                         path_model_cdf_test, output)
                    #            shrp.module1(path_emission_cdf_test, nc_redarea,
                    #                         path_reduction_txt, path_model_cdf_test,
                    #                         output)
                stop = time()
                print('Module 1 for {}, {}, run time: {} sec.'.format(country, prec,
                      (stop-start)))
                remove(proglog_filename)
            else:
                print('WARNING: previously calculated results will be used')

            # Get results
            deltaconc = output + 'delta_concentration.nc'
            rootgrp = Dataset(deltaconc, mode='r')
            d_pm25_conc = rootgrp.variables['delta_concentration'][:]
            rootgrp.close()
            path_areacountry_nc = 'workdir\\area_{}.nc'.format(country)
            rootgrp = Dataset(path_areacountry_nc, mode='r')
        #    areacountry = rootgrp.variables['AREA'][:]
            # bc_pm25_units = rootgrp.variables['conc'].units
            rootgrp.close()
            # Substitute 0s to NaN
            conca =  np.where(np.isnan(d_pm25_conc), 0, (d_pm25_conc))
#            print(df_red)
#            print(areatype, np.sum(conca))
            conc = conc + conca
        outfile=outputfin + 'delta_concentration.nc'
        if os.path.exists(outfile):
            print('removing file')
            os.remove(outfile) 
        write_nc(conc, outputfin + 'delta_concentration.nc', 'delta_concentration', 'microg/m3')

    rootgrp = Dataset(path_base_conc_cdf_test, mode='r')
    bc_pm25_conc = rootgrp.variables['conc'][:]
    # bc_pm25_units = rootgrp.variables['conc'].units
    rootgrp.close()
    bc_conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))
    # Total popluation in the country
    # poptot = np.sum(popall * areacountry / 100)
    exp = 1000 * conc * (popall) # ng p /m3
    bc_exp = 1000 * bc_conc * (popall)  # ng p /m3
    exp_per = (np.divide(exp, bc_exp, out=np.zeros_like(exp), where=bc_exp!=0))*100
    
    outfile = outputfin + 'exp{}{}{}{}{}per.nc'.format(country,name,order,an_ty,opt)
    if os.path.exists(outfile):
        print('removing file')
        os.remove(outfile)   
    write_nc(exp_per, outfile, 'exposure', '%')

    return exp_per

def save_xls(dic_dfs, xls_path):
    writer = ExcelWriter(xls_path)
    for n, df in enumerate(dic_dfs):
        dic_dfs[df].to_excel(writer,'%s' % df)
    writer.save()

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
    "text.usetex": True, #True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
#    'savefig.bbox': None,
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

def f1(x):
    return '%.f' % x

if __name__ == '__main__':

    # General parameters
    sources = ['BE', 'DE','ES', 'FR', 'IT', 'UK']
    yearred = '2030'  # or 2030
    year = 2005  # according to the new NEC directive
    opts = ['I', 'II']  # I to order sectors in terms of eff (I is used in the article)
                        # II to order sectors in terms of potential      
    an_tys = ['sec', 'spa']  # 'sec' or 'spa' anal. type sectorial or spatial
    ords = ['inc', 'dec'] # orders 'inc' or 'dec' increasing or decreasing
                          # values of indicators
                          
    name = 'exposure' # esposure or concentration
    path_results = 'spaflextest3\\'
    path_figures = path_results
    
    # select emission inventory
    em_inv = 'GAINS' # 'EMEP', 'TNO', 'GAINS'

    # Population
    path_tiff = 'input/pop/7km_Qi_2010.tif'
    popall = tiftogridgeneral(path_tiff)

    # Area of each cell in the domain
    path_surf_nc = 'input/JRC01.nc'
    rootgrp = Dataset(path_surf_nc, mode='r')
    surf = rootgrp.variables['surface'][:]
    area_units = rootgrp.variables['surface'].units
    rootgrp.close()
    
    # Longitude and latitutde arrays
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()
    n_lon = len(longitude_array)  
    n_lat = len(latitude_array)  

    
    precursors = ['PPM', 'NOx', 'SOx', 'NH3', 'NMVOC']
    
    # Absolute value of emission reduction by country and precursor
    df_e_red_eu = pd.DataFrame(columns=precursors)
    df_e_inv_eu = read_inventory(year, em_inv, precursors)

    # Dataframe with the emission reductins per country according to NEC
    df_nec = pd.read_csv(path_results + 'NECdirective.csv', index_col=[0],
                         header=[0, 1])
    
#------------------------------------------------------------------------------
    # Start calculations: 
        
    opt = 'I'
    for an_ty in an_tys:
        nec_app(an_ty, opt, name, sources, path_results, path_figures, em_inv)

    grid_per_ch = {}
    name = name
    for country in sources:
        grid_per_ch[country] = {}
        for an_ty in an_tys:
            grid_per_ch[country][an_ty] = {}
            grid_per_ch[country][an_ty]['dec'] = perexppercountry(name,an_ty, opt, country, 'dec')
            grid_per_ch[country][an_ty]['inc'] = perexppercountry(name,an_ty, opt, country, 'inc')

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
    # REVIEWED 07/12/2017
    
    # Cities analysis
    path_bc_conc_nc = path_base_conc_cdf_test
    rootgrp = Dataset(path_bc_conc_nc, mode='r')
    bc_pm25_conc = rootgrp.variables['conc'][:]
    # bc_pm25_units = rootgrp.variables['conc'].units
    rootgrp.close()
    
    conc = np.where(np.isnan(bc_pm25_conc), 0, (bc_pm25_conc))

    cities = pd.read_excel(path_results + '150fuas_lat-lon_eng.xls',
                           index_col='URAU_CODE')
    measurements = pd.read_table(path_results + 'AA_airbase\\airbase_ereporting_2015_stations_within_fua150.txt',sep=';', skipinitialspace=True)
    measurements.stationtype=measurements.stationtype.str.strip()
    measurements.pollutant=measurements.pollutant.str.strip()

    for uraucode in cities.index:
        lat = cities['lat'].loc[uraucode]
        lon = cities['lon'].loc[uraucode]
        en_name = cities['ENGLISH_NAME'].loc[uraucode]
            # convert latitude and longitude string in float
        lat = float(lat)
        lon = float(lon)
            # get row index of latitude and col index of longitude
        i_lat_target = 0
        lat_error = float('inf')
        for i in range(len(latitude_array)):
            lat_dist = abs(lat - latitude_array[i])
            if lat_dist < lat_error:
                lat_error = lat_dist
                i_lat_target = i
        
        i_lon_target = 0
        lon_error = float('inf')
        for i in range(len(longitude_array)):
            lon_dist = abs(lon - longitude_array[i])
            if lon_dist < lon_error:
                lon_error = lon_dist
                i_lon_target = i
                
#        print(conc[i_lat_target, i_lon_target])   
        cities.loc[uraucode, 'bc_conc'] = conc[i_lat_target, i_lon_target]
#        cities['bc_conc']
        mes_conc = measurements['aqvalue'].loc[(np.where((measurements.urau_code.values == uraucode) & 
                     (measurements.stationtype.values == 'Background') & 
                     (measurements.pollutant.values == 'Particulate matter < 2.5 µm (aerosol)')))].mean()
        cities.loc[uraucode, 'mes_conc'] = mes_conc  
        country=uraucode[0:2]
        if country in sources:
            for an_ty in an_tys:
                for order in ords:
                    output = path_results + 'NEC_{}{}{}{}{}\\'.format(country,name,order,an_ty,opt)
                    dc_path = output + 'delta_concentration.nc'
                    rootgrp = Dataset(dc_path, mode='r')
                    dc_pm25_conc = rootgrp.variables['delta_concentration'][:]
                    rootgrp.close()
                    dconc = np.where(np.isnan(dc_pm25_conc), 0, (dc_pm25_conc))
                    cities.loc[uraucode, (order+'_'+an_ty)] = dconc[i_lat_target, i_lon_target] 
                    cities.loc[uraucode, (order+'_'+an_ty+'scaled')] = dconc[i_lat_target, i_lon_target]/conc[i_lat_target, i_lon_target]*cities.loc[uraucode, 'mes_conc']

    # saving updated cities database:
    cities.to_csv(path_results + 'AA_airbase\\data.csv') 
    
    # plotting cities resutls
    cities_plot = ['Bruxelles', 'Antwerpen', 'Liege', 
                   'Berlin', 'Hamburg', 'Koln', 
                   'Madrid', 'Barcelona', 'Valencia',
                   'Paris', 'Lyon', 'Marseille', 
                   'Roma', 'Milano', 'Torino',
                   'London', 'Manchester', 'Liverpool']
    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()
#    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_ylim([0,30])
#    ax = plt.minorticks_on()
    ax.tick_params(axis='x',which='minor',bottom='off')
    ind = np.arange(len(cities_plot))
#    ax.set_xticklabels(list(cities['cityname'].values), rotation=90)
    plt.xticks(ind, cities_plot, rotation=90)
    who_x = ind
    who_y = np.ones(len(ind))* 10
    aqd_y = np.ones(len(ind))* 20 # to check
    plt.plot(who_x, who_y, color='g', linestyle = 'dashed', label = 'WHO limit')
    plt.plot(who_x, aqd_y, color='r', linestyle = 'solid', label = 'AQD limit')
    patches = []
    citiesni = cities.set_index(['ASCII_NAME'])
    for ind, x in enumerate(cities_plot):
        print(ind, citiesni['mes_conc'].loc[x])
        plt.plot(ind, citiesni['mes_conc'].loc[x],'o',color = 'r')
        plt.plot(ind, citiesni['mes_conc'].loc[x]-citiesni['dec_secscaled'].loc[x],'D', color = 'b')
        plt.plot(ind, citiesni['mes_conc'].loc[x]-citiesni['inc_secscaled'].loc[x],'D',  mfc='none', color ='b')
        plt.plot(ind, citiesni['mes_conc'].loc[x]-citiesni['dec_spascaled'].loc[x],'v', color = 'g')
        plt.plot(ind, citiesni['mes_conc'].loc[x]-citiesni['inc_spascaled'].loc[x],'v', mfc='none', color = 'g')
    plt.ylabel('Concentration [$\mu$ g/m$^3$]')
    marker = mlines.Line2D([], [], color='b', marker='D',
                          label='sec. dec.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='b', marker='D', mfc='none',
                          label='sec. inc.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='g', marker='v',
                          label='spa. dec.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='g', marker='v', mfc='none',
                          label='spa. inc.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='o',
                          label='base case', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='None',
                          label='AQD limit', linestyle = 'solid')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='None',
                          label='WHO limit', linestyle = 'dashed')
    patches.append(marker)
    plt.legend(handles=patches, ncol=2, loc='best')
    fig.savefig(path_figures + 'cities.png', bbox_inches = "tight", dpi=300)
    fig.savefig(path_figures + 'cities.pgf', bbox_inches = "tight")
    fig.savefig(path_figures + 'cities.pdf', bbox_inches = "tight")
    plt.show()
    
    ## test
    cities.loc['FR001L1']
    mes_conc = measurements['aqvalue'].loc[(np.where((measurements.urau_code.values == 'FR001L1') & 
                     (measurements.stationtype.values == 'Traffic') & 
                     (measurements.pollutant.values == 'Particulate matter < 2.5 µm (aerosol)')))]
 
#------------------------------------------------------------------------------     
    # REVIEWED 11/12/2017
    # tables with reducctions 
    # Dataframe with the emission reductins per country according to NEC
    df_nec = pd.read_csv(path_results + 'NECdirective.csv', index_col=[0],
                         header=[0, 1])
    # Dataframe with the emissions per country according to the inventory
    df_e_inv_eu = read_inventory(2005, em_inv, precursors)*1000
    redper={}
    red={}
    yearslong=[2005,2010]*len(precursors)
    preclong=np.repeat(precursors,2)
    redpertable = pd.DataFrame(index=sources, columns=zip(yearslong,preclong))
    for precursor in precursors:
        # Calculate the absolute value of emissions after implementation of NEC:
        df_e_red_eu[precursor] = (1 - df_nec[yearred][precursor].ix[sources]) * df_e_inv_eu[precursor].ix[sources]
        # Absolute value of emission reduction from year 2010
        emi2010 = (read_inventory(2010, em_inv, [precursor])*1000).loc[sources][precursor]
        red[precursor] = emi2010 - df_e_red_eu[precursor]
        redper[precursor] = red[precursor]/emi2010 * 100 
        redpertable[2005,precursor]=df_nec[yearred][precursor].ix[sources]*100 
        redpertable[2010,precursor]=redper[precursor]
        
    print(redpertable.to_latex(float_format=f1)) 

#------------------------------------------------------------------------------        
#------------------------------------------------------------------------------
   


    target_cell_lat = 40.40625
    target_cell_lon = -3.6875

    # convert latitude and longitude string in float
    target_cell_lat = float(target_cell_lat)
    target_cell_lon = float(target_cell_lon)
    

    
    country = 'ES'
    order = 'dec'
    an_ty = an_tys[0]
    output = path_results + 'NEC_{}{}{}{}{}\\'.format(country,name,order,an_ty,opt)
    dc_path = output + 'delta_concentration.nc'
    rootgrp = Dataset(dc_path, mode='r')
    dc_pm25_conc = rootgrp.variables['delta_concentration'][:]
    bc_pm25_units = rootgrp.variables['delta_concentration'].units
    rootgrp.close()
    
    dconc = np.where(np.isnan(dc_pm25_conc), 0, (dc_pm25_conc))
    ac = conc - dconc
    ac[i_lat_target, i_lon_target]

#    cities['NUTS0'] = [city[:2] for city in cities['FUA_CODE']]
#    citysel = cities.set_index('NUTS0').loc[sources]
#    av_conc_new = {}
#    av_conc = {}
#    for city in citysel['FUA_CODE']:
#        print(city)
#        av_conc_new[city] = {}
##        print(city)
    
    area_city = gridint_toarray('FUA_CODE', 'parea', 'ES001L2')
    area_city = gridint_toarray('GCITY_CODE', 'parea', 'ES001C1')
    av_conc = (np.sum((conc * area_city * surf / 100) / np.sum(area_city * surf / 100)))
    print(av_conc)
    
    for an_ty in an_tys:
        av_conc_new[city][an_ty] = {}
        for order in ords:
            output = path_results + 'NEC_{}{}{}{}{}\\'.format(country,name,order,an_ty,opt)
#    output = path_results + '\\NEC_{}{}{}{}{}\\'.format(country,name,'dec',an_ty,opt)
            dc_path = output + 'delta_concentration.nc'
            rootgrp = Dataset(dc_path, mode='r')
            dc_pm25_conc = rootgrp.variables['delta_concentration'][:]
# bc_pm25_units = rootgrp.variables['conc'].units
            rootgrp.close()
            dconc = np.where(np.isnan(dc_pm25_conc), 0, (dc_pm25_conc))
            ac = bc_pm25_conc - dconc
            av_conc_new[city][an_ty][order] = (np.sum((ac * area_city * surf / 100) /
            np.sum(area_city * surf / 100)))
            print(av_conc_new[city][an_ty][order])
    cities = pd.read_csv(
                   path_results + 'cities.csv', index_col=[0])
    plt.clf()
    fig= plt.figure(figsize=figsize(1))
    ax = fig.add_subplot(111)
    ax.minorticks_on()
#    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_ylim([0,30])
#    ax = plt.minorticks_on()
    ax.tick_params(axis='x',which='minor',bottom='off')
    ind = np.arange(len(cities['cityname'].values))
#    ax.set_xticklabels(list(cities['cityname'].values), rotation=90)
    plt.xticks(ind, cities['cityname'].values, rotation=90)
    who_x = ind
    who_y = np.ones(len(ind))* 10
    aqd_y = np.ones(len(ind))* 20 # to check
    plt.plot(who_x, who_y, color='r', linestyle = 'dashed', label = 'WHO limit')
    plt.plot(who_x, aqd_y, color='r', linestyle = 'solid', label = 'AQD limit')
    patches = []
    for ind, x in enumerate(cities.index):
        print(ind, av_conc[x])
        plt.plot(ind, av_conc[x],'o',color = 'r')
        plt.plot(ind, av_conc_new[x]['sec']['dec'],'D', color = 'b')
        plt.plot(ind, av_conc_new[x]['sec']['inc'], 'D',  mfc='none', color ='b')
        plt.plot(ind, av_conc_new[x]['spa']['dec'],'v', color = 'g')
        plt.plot(ind, av_conc_new[x]['spa']['inc'], 'v',  mfc='none', color ='g')
    plt.ylabel('Average concentration [$\mu$ g/m$^3$]')
    marker = mlines.Line2D([], [], color='b', marker='D',
                          label='sec. dec.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='b', marker='D', mfc='none',
                          label='sec. inc.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='g', marker='v',
                          label='spa. dec.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='g', marker='v', mfc='none',
                          label='spa. inc.', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='o',
                          label='base case', linestyle = 'None')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='None',
                          label='AQD limit', linestyle = 'solid')
    patches.append(marker)
    marker = mlines.Line2D([], [], color='r', marker='None',
                          label='WHO limit', linestyle = 'dashed')
    patches.append(marker)
    plt.legend(handles=patches, ncol=2, loc='best')
    fig.savefig(path_figures + 'cities.png', bbox_inches = "tight", dpi=300)
    fig.savefig(path_figures + 'cities.pgf', bbox_inches = "tight")
    fig.savefig(path_figures + 'cities.pdf', bbox_inches = "tight")
    plt.show()

#    A= np.sum((1000 * conc * areacountry * surf / 100) / np.sum(areacountry * surf / 100))
#    b= np.sum((1000 * bc_conc * areacountry * surf / 100) / np.sum(areacountry * surf / 100))
#    C= np.sum((1000 * conc * areacountry * popall / 100) / np.sum(areacountry * popall / 100))
#    D= np.sum((1000 * bc_conc * areacountry * popall / 100) / np.sum(areacountry * popall / 100))


    # Draw gridded impacts
    plt.close('all')
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()
    lon, lat = np.meshgrid(lon_array, lat_array)
#    dfmat=df2mat(dc_dic)
    fig= plt.figure(figsize=figsize(1))
    m = Basemap(resolution='i',projection='cyl',
            llcrnrlon=min(lon_array), llcrnrlat=min(lat_array),
            urcrnrlon=max(lon_array), urcrnrlat=max(lat_array))#,
    #        llcrnrlon=min(lons), llcrnrlat=min(lats),
    #        urcrnrlon=max(lons), urcrnrlat=max(lats),
#            area_thresh = 0.1,lat_0=lat_0,lon_0=lon_0)
    m.drawcoastlines(linewidth=1.25);
    #m.drawstates()
    m.drawcountries(linewidth=1.25);
    # draw filled contours.
    x, y = m(lon, lat)
    clevs = np.arange(5,55,5)#[5,10,15,20,25,30,35,40,45,50,55,60]
#    cs = m.contourf(x,y,dfmat,clevs);
    for an_ty in an_tys:
        for order in ords:
            grid_per_ch_tot = np.sum(grid_per_ch[country][an_ty][order] for country in ['ES', 'DE'])  # to check
            cs = m.contourf(x,y,(grid_per_ch_tot[0]), clevs, cmap='jet', vmax=50, extend='max')
            m.bluemarble()
            m.drawcoastlines(linewidth=0.5);
            m.drawcountries(linewidth=0.5);
            cbar = m.colorbar(cs,location='right', pad="5%", label = 'impact reduction [\%]')
            cbar.ax.tick_params(labelsize=7)
            plt.savefig(path_figures + '{}{}{}{}.pdf'.format(name, an_ty, opt, order), bbox_inches = "tight", format='pdf')
#    plt.savefig(path_figures + '{}{}{}{}.pgf'.format(name, an_ty, opt, order), format='pgf')
            plt.show()

    def f1(x):
        return '%.2f' % x

    def f2(x):
        return '%.2f' % x


    for precursor in precursors:
        print(precursor)
        print(redms2[precursor].to_latex(float_format=f1))

    for precursor in precursors:
        print(precursor)
        print(df_releff[precursor].to_latex(float_format=f1))
    for precursor in precursors:
        print(precursor)
        print(redms2[precursor].to_latex(float_format=f1))
    for precursor in precursors:
        print(precursor)
#        print(redper.to_latex(float_format=f1))

    precursors = ['PPM', 'NOx', 'SOx', 'NH3', 'NMVOC']  
    # Absolute value of emission reduction by country and precursor
    df_e_red_eu = pd.DataFrame(columns=precursors)
    df_e_inv_eu = read_inventory(year, em_inv, precursors)
    df_e_inv_eu10 = read_inventory(2010, em_inv, precursors)
    # Dataframe with the emission reductins per country according to NEC
    df_nec = pd.read_csv(path_results + 'NECdirective.csv', index_col=[0],
                         header=[0, 1])
    
    print(df_e_inv_eu.loc[sources].to_latex(float_format=f1))
    print(df_e_inv_eu10.loc[sources].to_latex(float_format=f1))
#  DRAW FUAS
    rootgrp = Dataset(path_model_cdf_test, 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    rootgrp.close()
    lon, lat = np.meshgrid(lon_array, lat_array)
#    dfmat=df2mat(dc_dic)
    plt.close('all')
    fig= plt.figure(figsize=figsize(1))
    m = Basemap(resolution='i',projection='cyl',
            llcrnrlon=min(lon_array), llcrnrlat=min(lat_array),
            urcrnrlon=max(lon_array), urcrnrlat=max(lat_array))#,
    #        llcrnrlon=min(lons), llcrnrlat=min(lats),
    #        urcrnrlon=max(lons), urcrnrlat=max(lats),
#            area_thresh = 0.1,lat_0=lat_0,lon_0=lon_0)
    m.drawcoastlines(linewidth=0.5);
    m.drawcountries(linewidth=0.5);
    x, y = m(lon, lat)
    clevs = np.arange(5,105,5)#[5,10,15,20,25,30,35,40,45,50,55,60]
    path_fua_country = 'D:\\sherpa.git\\Sherpa\\asstest\\workdirreg\\'
    area = np.zeros((1, len(lat_array), len(lon_array)))
    for country in sources:
        rootgrp = Dataset(path_fua_country + 'area_fuas_{}.nc'.format(country), 'r')
        area_co=rootgrp.variables['AREA'][:]
        area = area + area_co
    cs = m.contourf(x,y,(area[0]), clevs)#, cmap='jet', vmax=50, extend='max')
    m.shadedrelief()
#    m.bluemarble(ax=None, scale=1)
    cbar = m.colorbar(cs,location='right', pad="5%", label = 'cell area [\%]')
    cbar.ax.tick_params(labelsize=7)
    plt.savefig(path_figures + 'fuas3.pdf', bbox_inches = "tight", format='pdf')
    plt.show()


