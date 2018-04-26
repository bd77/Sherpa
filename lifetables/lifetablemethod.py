# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 15:22:19 2018

Script to calculate the health impact assessment using life tables - 
it calculates the loss of life expectancy of a cohort over its lifetime. 
The script is flexible enough to allow (with little manual modifications)
the analysis of a scenario of concentration reductions over the lifetime of 
the cohort. 

This script allows to obtain the same order of magnitude values as in TSAP 11


input:
    - country
    - l_int_age: minimum age interval of cohort to be consider in the analysis
      can be chosen within this list (example consider the cohort that today 
      is '< 1 year' and older):       
       '< 1 year', '1 year', '2 years', '3 years', '4 years', '5 - 9',
       '10 - 14', '15 - 19', '20 - 24', '25 - 29', '30 - 34', '35 - 39',
       '40 - 44', '45 - 49', '50 - 54', '55 - 59', '60 - 64', '65 - 69',
       '70 - 74', '75 - 79', '80 - 84', '85 - 89', '90 - 95', '95 +'
    - age_impact: minimum age on which the impact is assessed (e.g. 30 years
      and older)
        0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65
        70, 75, 80, 85
    - path_mortbaseline - path to file with mortality baseline (see below)
    - path_population - path to file with population baseline (see below)
            
prepare data: 
        # go to: 
        http://apps.who.int/healthinfo/statistics/mortality/causeofdeath_query/
        # for mortality: 
        Extract data by country, year, sex and age for a few selected causes of death by coding systems
            # select: Detailed ICD-
            # select for group A the causes that you will modify 
            # select for group B the total causes
            # select all age groups
            # save data as csv - NEVER OPEN THIS WITH EXCEL (or you will get dates)
        # for population data
        Extract population data by country, year, sex and age
        
References:
[1] B.G. Miller and J.F. Hurley, Life table methods for quantitative impact 
    assessment in chronic mortality (2002)

[2] World Health Organization Europe, “Quantification of the Health Effects of
     Exposure to Air Pollution,” 2001.
     
[3] World Health Organization Europe, 2013. Health risks of air pollution
    in Europe - HRAPIE project - Recommendations for concentration–response
    functions for cost–benefit analysis of particulate matter, ozone and
    nitrogen dioxide, Copenhagen Ø, Denmark.

[4] http://apps.who.int/healthinfo/statistics/mortality/causeofdeath_query/

[5] World Health Organization Europe, 2017. AirQ+: software tool for health
    risk assessment of air pollution.
    
[6] Description of the Eurostat method for the calculation of the life expectancies
    at all ages. European Commission, Eurostat - Directorate F, Unit F-2 population

@author: peduzem
"""

import pandas as pd
import numpy as np
import sys

def lifetable(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age):    
    '''Function that produces calculates over time the life years and the 
        deaths of the population alive in the initial cohort calculating the 
        equivalent of 'table 5' in [1]
       input:
           - df_table_h_mod: dataframe that modifies the baseline health 
            impact rates (e.g. all cause mortality) in terms of a fraction
            (e.g. 0.99 corresponds to a 1% reduction of the baseline
            risk rate).
           - h_table: baseline of health impact of interest
           - h_rest_table: rest of health impact (TOTAL cause mortality)
           - e: entry population 
           - country
           - keys: keys of age ranges 
           - values: corresponding age values 
           - dates: corresponding calendar dates
           - l_int_age: minimum age interval of cohort to be consider in the 
           analysis
        output: 
            - life_years_millions: million of life years lived by the cohort
            - klife_years_per100k: thousands of live years lived per 100kppl
    '''
    map_age = dict((zip(keys, values)))
    
    # corresponds to table 4 in [1] risk rates over their lifetimes
    # HP: baseline riskrates in future dates are the same as the initial 
    # cohorts ones dipending on age. 
    df_table4h=pd.DataFrame(columns=dates, index=values)
    for date in dates:
        df_table4h.loc[:,date]=h_table[keys].loc[country].values


    df_table4h_rest=pd.DataFrame(columns=dates, index=values)
    for date in dates:
        df_table4h_rest.loc[:,date]=h_rest_table[keys].loc[country].values

    # modified risk rates
    df_table4h_ac =  df_table4h*df_table_h_mod + df_table4h_rest 
    
    # length of age groups
    dvalue = [values[i+1]-values[i] for i in range(len(values)-1)] + [5]
    
    # survival probability over dvalue years
    df_table4s_ac= ((2-df_table4h_ac)/(2+df_table4h_ac)).pow(dvalue, axis='rows') 
    
    # impose survivel 0 for last age
    df_table4s_ac.loc[values[-1]]=0          
    
    # transform dataframe as matrix
    m_s_ac = df_table4s_ac.as_matrix() 
    
    # initialise matrices for results (corresponding to table 5)
    m_d=np.zeros(np.shape(m_s_ac)) # deaths
    m_L=m_d # life years lived
    
    for i, age in enumerate(values): 
        # get each diagonal of each age group and compute the life tables:
        h_diag = np.diag(df_table4h_ac.as_matrix(), k=-i)
        s_diag = np.diag(df_table4s_ac.as_matrix(), k=-i)
        h_diag = h_diag.copy()

        # change death rate of last group to 1
        h_diag[-1]=1
    
        # cumulative survival probability 
        cum_s_diag = np.cumprod(s_diag)
        
        # survivors at age (for each year)
        l_diag =  [1] + [cum_s_diag[it] for it in range(len(cum_s_diag)-1)] 

        # people of the inital cohort surving for each future age group
        l_diag = list(np.array(l_diag)*e[age].loc[country])
        
        # deaths occuring in each age group 
        d_diag = np.array(l_diag)-np.array(l_diag[1:]+[0])  
        
        # years lived in the future for each 
        L_diag = list(np.array(list((np.array(l_diag[1:]) + 0.5 * d_diag[0:-1])))*dvalue[i:-1]) + [l_diag[-1]/h_diag[-1]]
        
        m_d= m_d + np.diag(d_diag, -i)
        m_L= m_L + np.diag(L_diag, -i)

    ml=[m_d, m_L]

    min_age=map_age[l_int_age]
    ind = values.index(min_age)
    
    # results dictionary 
    r={}
    results_keys=['total deaths', 'total years of life lived']

    for i,m in enumerate(ml):

        rows=np.shape(m[ind:,:])[0]
        r[results_keys[i]]=np.array([np.trace(m[ind:,:], x) 
                       for x in range(-(rows-1), 1)]).sum()
   
    checksum=r[results_keys[0]] /e[values[ind:]].loc[country].sum()   
    if checksum == 1:
       print('ALL GOOD: all initial people eventually die')
    else: 
       print('WARNING: there is a problem, ratio between deaths and initial people is: ', checksum)

    life_years_millions = r[results_keys[1]]/(10**6)
    klife_years_per100k = (100*10**3)*(r[results_keys[1]]/10**3)/(e[values[ind:]].loc[country].sum())   
    return life_years_millions, klife_years_per100k



if __name__ == '__main__':
   
    # input
    country = 'Italy'
    l_int_age = '< 1 year' # minimum age of cohort to consider in the analysis '25 - 29'#
    age_impact = 30 # minimum age on which the health impact is applied
    avg_exp = 8.359

    # average exposure change 
    
    # input files paths
    path_mortbaseline = 'D://sherpa.git//Sherpa//lifetables//extract2allcountries.csv' 
    path_population = 'D://sherpa.git//Sherpa//lifetables//populationallcountries.csv' 
    
    ### concentration response function
    # CONCENTRATION RESPONSE FUNCTION:
    # From [2] Table 1 
    # Estimate of mortality, all-cause (natural) age 30+ years
    # PM2.5 Annual mean 
    # RR = 1.062 (1.04-1.083) 95% CI per 10 microg/m3 

    # From [5] Beta: 
    lbeta = 0.003922071315328133  # lower bound
    mbeta = 0.006015392281974714  # average value
    hbeta = 0.007973496801885352  # higher bound
    # RR = e^(beta x)=(e^(beta*10))^(x/10) = 1.062^(x/10) CI = 1.04, 1.083
    # (we obtain the same values reported in [2])


    # exponential approximation
    # AF = (RR -1)/RR = e^bx -1 / (e^bx) = 1 -e^(-bx) 
    beta = [lbeta, mbeta, hbeta]
    pt = len(beta)
    af =  np.array([(1-(np.exp(-beta[i]*avg_exp))) for i in range(len(beta))])

    
    ### This part should be adapted to read the data 
    df_mort = pd.read_csv(path_mortbaseline, skiprows=40, index_col=[0,3])
    frmt_mort=list(set(df_mort['Format'].loc[country].loc['B'].values))[0]

    year_mort=list(set(df_mort['Year'].loc[country].loc['B'].values))[0]
    df_mort.drop(['ICD rev and list','Sex', 'Unknown', 'Sum of selected ages','Format', 'Year'], axis=1, inplace=True)

    # dataframe with all cause mortality (excluding external causes)     
    df_mort_A = df_mort.xs('A', level=1).sum(level=0)

    # dataframe with total cause mortality (including external causes)         
    df_mort_B=df_mort.xs('B', level=1).sum(level=0)
    
    # dataframe with mortality that are not impacted by the cause(s) of choice         
    df_mort_BA=df_mort_B-df_mort_A
    
    # population data    
    df_pop = pd.read_csv(path_population, skiprows=38, index_col=[0])
    frmt_pop=list(set(df_pop['Format'].loc[country].values))[0]
    
        
    year_pop=list(set(df_pop['Year'].loc[country].values))[0]
    
    df_pop = df_pop.sum(level=0) # summing F and M
    df_pop = df_pop[df_mort.columns] # taking only columns of interest
    
    if year_pop == year_mort: 
        print('OK: population and mortality data are both from year: ', year_pop)
    else: 
        print('WARNING: inconsistent population and mortality years')                  
    year=year_pop

    if frmt_mort==0 and frmt_pop==0:
    # time frames (years, dates, etc.)
        keys = list(df_pop.columns)
        values = [0, 1, 2, 3, 4] + list(np.arange(5,100,5))
    elif frmt_pop==1 or frmt_mort==1: 
        keys = list(df_pop.columns)[:-2]
        values = [0, 1, 2, 3, 4] + list(np.arange(5,90,5))
        df_pop.loc[country,'85 - 89'] = df_pop[['85 - 89','90 - 95', '95 +']].loc[country].sum()
        df_mort_A.loc[country,'85 - 89'] = df_mort_A[['85 - 89','90 - 95', '95 +']].loc[country].sum()
        df_mort_BA.loc[country,'85 - 89'] = df_mort_BA[['85 - 89','90 - 95', '95 +']].loc[country].sum()
        df_mort_B.loc[country,'85 - 89'] = df_mort_B[['85 - 89','90 - 95', '95 +']].loc[country].sum()
    else:
        sys.exit('NOT YET IMPLEMENTED FOR THE FORMAT')

    dates= [value+year for value in values]

    # entry population
    e = df_pop[keys] + df_mort_B[keys]/2
    e.columns=values     
    
    # calculating mortality rate          
    h_table=pd.DataFrame(columns=keys) # mortality rate for each cause of interest
    h_rest_table=pd.DataFrame(columns=keys) # rest of mortality rate
    h_table[keys]= df_mort_A[keys]/df_pop[keys] 
    h_rest_table[keys]= df_mort_BA[keys]/df_pop[keys]
    
    # table that modifies the risk for the causes of choice
    df_table_h_mod=pd.DataFrame(1, columns=dates, index=values)
    bl_life_years_millions, bl_klife_years_per100k = lifetable(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age)
    
    df_table_h_mod=pd.DataFrame(1, columns=dates, index=values)

    df_table_h_mod=pd.DataFrame(1, columns=dates, index=values)
    
    # example: 
#    df_table_h_mod.loc['30':, 2063:]=0.85
#    df_table_h_mod.loc[:]=0.99

    df_table_h_mod.loc[age_impact:]=1-af[1]
    df_table_h_mod.to_excel('D://sherpa.git//Sherpa//lifetables//h_mod.xls')
    life_years_millions, klife_years_per100k = lifetable(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age)
    
    print('Total life years gained/lost (millions)',life_years_millions-bl_life_years_millions )
    print('Life years gained (thousands)\n',
          'per 100 000 population:', klife_years_per100k-bl_klife_years_per100k )
    print('Months gained per person:', ((klife_years_per100k-bl_klife_years_per100k)*1000)*12/100000 )