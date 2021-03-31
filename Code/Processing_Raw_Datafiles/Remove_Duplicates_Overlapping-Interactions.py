# -*- coding: utf-8 -*-
"""
This code removes the duplicates in the datafile 'Overlap_First_Interactors.xlsx',
such that each gene name occurs only in the colum 'gene' OR the column 'Interactor'
(to avoid double counting) and it assigns a single interaction type between two 
genes (Positive, Negative or Dubious)

This calculation requires the following data file: 
    Overlap_First_Interactors.xlsx
    
This code outputs the following data file:
    Overlap_First_Interactors_NoDuplicates.xlsx
"""

import pandas as pd

Interaction_Matrix = pd.read_excel (r'Input_Data_Files\Overlap_First_Interactors.xlsx')

Genes = pd.unique(Interaction_Matrix['Gene'])
Interaction_Matrix_New = pd.DataFrame(data=None, columns=Interaction_Matrix.columns)

for i in Genes:
    
    Query = Interaction_Matrix.loc[Interaction_Matrix['Gene'] == i,:]
    Interactors = pd.unique(Query['Interactor'])
    
    for j in Interactors:
        
        Int_Types = Query.loc[Query['Interactor']==j,'Interaction_Type']
        
        if (sum(Int_Types == 'Negative Genetic') == Int_Types.size) | (sum(Int_Types == 'Positive Genetic') == Int_Types.size): 
            
            Interaction_Matrix_New = Interaction_Matrix_New.append(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) \
                                                        & (Interaction_Matrix['Interactor'] == j)].iloc[0],ignore_index=True)
            
        else:
            
            Interaction_Matrix_New = Interaction_Matrix_New.append(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) \
                                                        & (Interaction_Matrix['Interactor'] == j)].iloc[0],ignore_index=True)
            
            Interaction_Matrix_New.loc[(Interaction_Matrix_New['Gene'] == i) & (Interaction_Matrix_New['Interactor'] == j),\
                                   'Interaction_Type'] = 'Dubious'
            

Interaction_Matrix_New.to_excel(r'Output_Data_Files\Overlap_First_Interactors_NoDuplicates.xlsx', index = False)