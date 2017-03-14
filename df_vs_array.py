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
