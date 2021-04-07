# -*- coding: utf-8 -*-
"""
This code determines the diversity of each interactor 
(how often it interacts with polarity vs diauxic shft genes)
This calculation requires the following data files: 
    Overlap_First_Interactors_NoDuplicates.xlsx

"""


import pandas as pd
import numpy as np
import seaborn as sns
import scipy
from scipy.stats import norm
#%% Import interaction matrix (in which each interactor occurs only once)
Interaction_Matrix = pd.read_excel (r'Data_Files/Overlap_First_Interactors_NoDuplicates.xlsx')

#%% Preparing datasets

Genes = pd.unique(Interaction_Matrix['Gene'])
Col_Vars = ['Gene','Number_of_Positive_Interactions_DiauxicShift','Number_of_Negative_Interactions_DiauxicShift'\
            ,'Number_of_Positive_Interactions_Polarity','Number_of_Negative_Interactions_Polarity']
Gene_Positive_Negative_Number = pd.DataFrame(data=None, columns=Col_Vars)

Interaction_Matrix.index=Interaction_Matrix['Gene']
Interaction_Matrix.drop(columns='Gene',inplace=True)


#%% # Loop over each gene to get the number of interactions it has with genes in 
# the polarity and genes in the diauxic shift modules
for i in Genes: 
    
    Numb_Pos_Int_DS = len(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) & \
                              (Interaction_Matrix['Interaction_Type'] == 'Positive Genetic') \
                              & (Interaction_Matrix['Interactor_Module'] == 'Diauxic Shift')])
                       
    Numb_Neg_Int_DS = len(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) & \
                              (Interaction_Matrix['Interaction_Type'] == 'Negative Genetic') \
                              & (Interaction_Matrix['Interactor_Module'] == 'Diauxic Shift')])
    
    Numb_Pos_Int_Pol = len(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) & \
                              (Interaction_Matrix['Interaction_Type'] == 'Positive Genetic') \
                              & (Interaction_Matrix['Interactor_Module'] == 'Polarity')])
                       
    Numb_Neg_Int_Pol = len(Interaction_Matrix[(Interaction_Matrix['Gene'] == i) & \
                              (Interaction_Matrix['Interaction_Type'] == 'Negative Genetic') \
                              & (Interaction_Matrix['Interactor_Module'] == 'Polarity')])
    
    Gene_Positive_Negative_Number = Gene_Positive_Negative_Number.append({'Gene': i, \
                                                                          'Number_of_Positive_Interactions_DiauxicShift':Numb_Pos_Int_DS,\
                                                                          'Number_of_Negative_Interactions_DiauxicShift':Numb_Neg_Int_DS,\
                                                                          'Number_of_Positive_Interactions_Polarity':Numb_Pos_Int_Pol,\
                                                                          'Number_of_Negative_Interactions_Polarity':Numb_Neg_Int_Pol},\
                                                                           ignore_index=True)
            


#%% Determine the distribution of the polarity/diauxic shift interactions for the gene collection

Gene_Positive_Negative_Number['Total_Diauxic'] = \
                                Gene_Positive_Negative_Number['Number_of_Negative_Interactions_DiauxicShift']+\
                                Gene_Positive_Negative_Number['Number_of_Positive_Interactions_DiauxicShift']

Gene_Positive_Negative_Number['Total_Polarity'] = \
                                Gene_Positive_Negative_Number['Number_of_Negative_Interactions_Polarity']+\
                                Gene_Positive_Negative_Number['Number_of_Positive_Interactions_Polarity']

Gene_Positive_Negative_Number['Polarity_over_Diauxic'] = (Gene_Positive_Negative_Number['Total_Polarity']/123-\
                                                        Gene_Positive_Negative_Number['Total_Diauxic']/38)\
                                                        /(Gene_Positive_Negative_Number['Total_Diauxic']/38+\
                                                         Gene_Positive_Negative_Number['Total_Polarity']/123)
# questions: 
    # - where do 123 and 38 come from? 
    # - can you explain this expresion you used to compute the ['polarity_over_diauxic']? 
#%% 
# can you explain the rationale behind these magnitudes? 

Dist_Overlap_Inter = Gene_Positive_Negative_Number[(Gene_Positive_Negative_Number['Total_Polarity']!=0) \
                                                  & (Gene_Positive_Negative_Number['Total_Diauxic']!=0)] 
Dist_Non_Overlap_Inter = Gene_Positive_Negative_Number[(Gene_Positive_Negative_Number['Total_Polarity']==0) \
                                                  | (Gene_Positive_Negative_Number['Total_Diauxic']==0)]
#%%
# Create a histogram of the 'Diversity'

sns.set(rc={'figure.figsize':(2.4,3)})
sns.set_style('whitegrid')
spacing = 0.05
n_bins = int(2//(2*spacing))

plot = sns.distplot(Dist_Overlap_Inter.Polarity_over_Diauxic.tolist(), kde=False, \
                    fit=scipy.stats.norm, norm_hist = True, bins = n_bins, vertical=True)

plot.set_xlabel("# of genes (scaled)",fontsize=14)
plot.set_ylabel("Diversity",fontsize=14) 
plot.tick_params(labelsize=11)

params = scipy.stats.norm.fit(Dist_Overlap_Inter.Polarity_over_Diauxic)
print(params)