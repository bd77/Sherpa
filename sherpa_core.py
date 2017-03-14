import numpy as np
import pandas as pd

'''
NAME
    Core of sherpa model calculation using dataframes and flatweight as an option
PURPOSE
    Calculate the delta conc for one precursor on all snaps, the output is a pd.Series which will be concatenated in a Df in the maim 
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
def sherpa_df(pr,model_idx,alldists,emissiondelta,inner_radius=False):
    emidelta=emissiondelta.loc[pr]
    #fastest
    window_all=(1.+alldists)**(-model_idx.loc['omega'][pr])
    dc_cells=model_idx.loc['alpha'][pr]*(emidelta*window_all) 
    if inner_radius is not False : #if not zero then use info on flatweight
        outer_cells=alldists[(alldists>inner_radius) & (emidelta.sum()>0)].index
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
    Core of sherpa model calculation using np array on 3D emissions without flat weight (for testing)
PURPOSE
    Calculate the delta conc for all precursors and snaps, output is a 3d array (snaps,emissing_points,precursors) 
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
def sherpa_np(model_idx,alldists_dic,emi_array):
    preclist=emi_array['dimnames']['precursors']

    #alldists_1=1.+alldists
    #window_all=np.power.outer(alldists_1, -np.array(model_idx.loc['omega'][preclist]))
    #window_all=np.asarray(list(map(lambda p: alldists_1**(-np.array(model_idx.loc['omega'][p])),preclist))).T
    #fastest
    alldists_1=1.+alldists_dic     
    window_all=np.array(pd.concat(list(map(lambda p: alldists_1**(-np.array(model_idx.loc['omega'][p])),preclist)),axis=1))

    dc_cells=np.multiply(np.multiply(emi_array['data'],window_all),np.array(model_idx.loc['alpha'][preclist]))
    return dc_cells  
    
'''
NAME
    Core of sherpa model calculation using np array and numexpr on 3D emissions without flat weight (for testing, fastest)
PURPOSE
    For one receptor grid point calculate the gridded delta conc from all emitting grid points for all precursors and all snaps without flatweight
    mdel_idx are the receptors SR coefficients
    alldists is a vector containing the distances of all emitting grids to the receptor 
    emi_array is a 3D array (snaps, emissing grid points, precursors)
    The output dc_cells is a 3D matrix (snaps,emitting grid points, precursors)
PROGRAMMER(S)
    Denise Pernigotti, March 2017
REVISION HISTORY
    
REFERENCES
    
'''
import numexpr as ne

def sherpa_ne(model_idx,alldists,emi_array):
    preclist=emi_array['dimnames']['precursors']
    emi=emi_array['data']
    window = np.zeros((len(alldists),len(preclist)))
    alpha=model_idx['alpha'][preclist]
    alldists_1=alldists+1
    for ip,prec in enumerate(preclist):
        omega_neg=-model_idx.loc['omega'][prec]
        window[:,ip]=ne.evaluate("alldists_1**omega_neg")
    dc_cells=ne.evaluate("(emi*window)*alpha")

    return dc_cells  

