# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 15:22:19 2018

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
        
References
B.G. Miller and J.F. Hurley, Life table methods for quantitative impact 
assessment in chronic mortality (2002)

http://data.euro.who.int/dmdb/ [Accessed December 13, 2016].
http://apps.who.int/healthinfo/statistics/mortality/causeofdeath_query/

@author: peduzem
"""

import pandas as pd
import numpy as np


def table5(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age):    
    
    map_age = dict((zip(keys, values)))
      
    df_table4h=pd.DataFrame(columns=dates, index=values)
    for date in dates:
        df_table4h.loc[:,date]=h_table[keys].loc[country].values


    df_table4h_rest=pd.DataFrame(columns=dates, index=values)
    for date in dates:
        df_table4h_rest.loc[:,date]=h_rest_table[keys].loc['Belgium'].values

    df_table4h_ac =  df_table4h*df_table_h_mod + df_table4h_rest 
    
    # length of age groups
    dvalue = [values[i+1]-values[i] for i in range(len(values)-1)] + [5]
    
    # survival
    df_table4s_ac= ((2-df_table4h_ac)/(2+df_table4h_ac)).pow(dvalue, axis='rows') 
    # impose survivel 0 for last age
    df_table4s_ac.loc[values[-1]]=0          
    
    # transform dataframe as matrix
    m_s_ac = df_table4s_ac.as_matrix() 
    
    # initialise matrices for results (corresponding to table 5)
    m_d=np.zeros(np.shape(m_s_ac))
    m_L=m_d
    
    for i, age in enumerate(values): 
        # get each diagonal of each age group and compute the life tables: 
#        i = 0
#        age = values[i]
        
        h_diag = np.diag(df_table4h_ac.as_matrix(), k=-i)
        s_diag = np.diag(df_table4s_ac.as_matrix(), k=-i)
        h_diag = h_diag.copy()
        # change death rate of last group to 1
        h_diag[-1]=1
    
        # cumulative survival probability 
        cum_s_diag = np.cumprod(s_diag)
        l_diag =  [1] + [cum_s_diag[i] for i in range(len(cum_s_diag)-1)] 
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
#    r_spec={}
#    results_keys_spec=['check deaths', 'total years of life lived']
    for i,m in enumerate(ml):
        rows=np.shape(m[ind:,:])[0]
        r[results_keys[i]]=np.array([np.trace(m[ind:,:], x) 
                       for x in range(-(rows-1), ind+1)]).sum()
   
    checksum=r[results_keys[0]] /e[values[ind:]].loc[country].sum()   
    if checksum == 1:
       print('ALL GOOD: all initial people eventually die')
    else: 
       print('WARNING: there is a problem, ratio between deaths and initial people is: ', checksum)

    life_years_millions = r[results_keys[1]]/(10**6)
    klife_years_per100k = (100*10**3)*(r[results_keys[1]]/10**3)/(e[values[ind:]].loc[country].sum())   
    return life_years_millions, klife_years_per100k



if __name__ == '__main__':
    
    ### 
    path_mortbaseline = 'D://sherpa.git//Sherpa//lifetables//extract2tot.csv' 
#    path_mortbaseline = 'D://sherpa.git//Sherpa//lifetables//extract2ew.csv' 
    # input
    country='Belgium'
    l_int_age='< 1 year' # minimum age of cohort to consider in the analysis
    
    ### This part should be adapted to read the data now working for format population 00
    
    df_mort = pd.read_csv(path_mortbaseline, skiprows=9, index_col=[0,3])
    df_mort.drop(['Format','ICD rev and list','Sex', 'Unknown'], axis=1, inplace=True)
    
    # dataframe with all cause mortality (or cause of choice)
    df_mort_A = df_mort.xs('A', level=1).sum(level=0)
    df_mort_A['Year']= df_mort_A['Year']/2
    
    # dataframe with total cause mortality (including external causes)         
    df_mort_B=df_mort.xs('B', level=1).sum(level=0)
    df_mort_B['Year']= df_mort_B['Year']/2
    
    # dataframe with mortality that are not impacted by the cause(s) of choice         
    df_mort_BA=df_mort_B-df_mort_A
    df_mort_BA['Year']=df_mort_B['Year']     
    
    # populatin data 
    path_population = 'D://sherpa.git//Sherpa//lifetables//population.csv' 
    df_pop = pd.read_csv(path_population, skiprows=6, index_col=[0])
    df_pop = df_pop.sum(level=0) # summing F and M
    df_pop['Year']= df_pop['Year']/2
    df_pop = df_pop[df_mort.columns] # taking only columns of interest
    
    year=df_pop['Year'].loc[country]
    keys = list(df_pop.columns[2:])
    values = [0, 1, 2, 3, 4] + list(np.arange(5,100,5))
    dates= [value+year for value in values]
    # entry population (to check**)
    e = df_pop[keys] + df_mort_B[keys]/2
    e.columns=values     
    
    # calculating mortality rate          
    h_table=pd.DataFrame(columns=keys) # mortality rate for each cause of interest
    h_rest_table=pd.DataFrame(columns=keys) # rest of mortality rate
    h_table[keys]= df_mort_A[keys]/df_pop[keys] 
    h_rest_table[keys]= df_mort_BA[keys]/df_pop[keys]
    
    # table that modifies the risk for the causes of choice
    df_table_h_mod=pd.DataFrame(1, columns=dates, index=values)
    bl_life_years_millions, bl_klife_years_per100k = table5(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age)
    df_table_h_mod=pd.DataFrame(1, columns=dates, index=values)
#    df_table_h_mod=pd.DataFrame(0.85, columns=dates, index=values)
    # example: 
    df_table_h_mod.loc['30':, 2063:]=0.85
#    df_table_h_mod.loc['30':]=0.85
    life_years_millions, klife_years_per100k = table5(df_table_h_mod, h_table, h_rest_table, e, country, keys, values, dates, l_int_age)

    print('Total life years gained (millions)',life_years_millions-bl_life_years_millions )
    print('Life years gained (thousands)\n',
          '100 000 population', klife_years_per100k-bl_klife_years_per100k )
