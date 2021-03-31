# -*- coding: utf-8 -*-
"""
This script calculates the percentage overlap between the first (neighbouring) 
genetic interactors of the polarity and diauxic shift modules. 

This calculation requires the following data files: 
    
    Genes_Diauxic_Shift_Description_Interactors.xlsx
    Genes_Diauxic_Shift_Description_Interactors_Polarity_Description.xlsx
    Polarity_Genes_Description_Interactors.xlsx
    Polarity_Genes_Description_Interactors_Diauxic_Shift_Description.xlsx
    Overlap_First_Interactors.xlsx
    Positive_Negative_Counts_PerGene.xlsx
"""
import pandas as pd
import matplotlib.pyplot as plt
import os

# Read in all the different data files and store them as dataframes

## Data related to the diauxic shift first interactors
DS_Interaction_Matrix = pd.read_excel(r'Data_Files\Genes_Diauxic_Shift_Description_Interactors.xlsx')
DS_Interaction_wPolarity_Matrix = pd.read_excel(r'Data_Files\Genes_Diauxic_Shift_Description_Interactors_Polarity_Description.xlsx')

## Data related to the polarity first interactors
Pol_Interaction_Matrix = pd.read_excel(r'Data_Files\Polarity_Genes_Description_Interactors.xlsx')
Pol_Interaction_wDS_Matrix = pd.read_excel(r'Data_Files\Polarity_Genes_Description_Interactors_Diauxic_Shift_Description.xlsx')

## Data related to the overlappig set of first interactors from the polarity 
## and diauxic shift interactors
Interaction_Pol_DS__Matrix = pd.read_excel(r'Data_Files\Overlap_First_Interactors.xlsx')


# Create dataframe containig all first interactors (polarity and diauxic shift combined)
All_Interactors = DS_Interaction_Matrix['Gene > Interactions > Participant 2 . Secondary Identifier']
All_Interactors = All_Interactors.append(Pol_Interaction_Matrix['Gene > Interactions > Participant 2 . Secondary Identifier'],ignore_index=True)

# Create dataframe containing all direct interactors (polarity gene <-> Diauxic shift gene)
Direct_Interactors = Pol_Interaction_wDS_Matrix['Gene > Systematic Name']
Direct_Interactors = Direct_Interactors.append(DS_Interaction_wPolarity_Matrix['Gene > Systematic Name'])

# Create dataframe containg all indirect interators 
# (Polarity gene <-> first interactor <-> Diauxic shift gene)
Indirect_Interactors = Interaction_Pol_DS__Matrix['Gene']

# Determine number of Diauxic shift (DS) and Polarity (Pol) interactors
Number_DS_Interact = len(pd.unique(DS_Interaction_Matrix['Gene > Interactions > Participant 2 . Secondary Identifier']))
Number_Pol_Interact = len(pd.unique(Pol_Interaction_Matrix['Gene > Interactions > Participant 2 . Secondary Identifier']))

# Determine total number of interactions, direct interactions and 
# indirect interactions between the polarity and diauxic shift modules
Total_Number_Interactors = len(pd.unique(All_Interactors))
Number_Direct_Interact = len(pd.unique(Direct_Interactors))
Number_Indirect_Interact = len(pd.unique(Indirect_Interactors))


# Create pie shart that shows the percentage of interactions that are 
# not shared/shared(= Direct interactions + Indirect interactions)
Not_Shared_Perc = (Total_Number_Interactors-Number_Indirect_Interact-Number_Direct_Interact)/Total_Number_Interactors*100
Direct_Perc = Number_Direct_Interact/Total_Number_Interactors*100
Indirect_Perc = Number_Indirect_Interact/Total_Number_Interactors*100

labels = 'Shared','Not Shared' 
sizes = [Direct_Perc+Indirect_Perc,Not_Shared_Perc]
cs1 = 'navajowhite','tomato','lightsteelblue'
explode = (0, 0.15)  # only "explode" the 2nd slice
fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', textprops={'fontsize': 26},
        shadow=False, startangle=90, colors = cs1)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.tight_layout()

# Create pie chart that shows the percentage of positive vs negative interactions 
# of those that are shared

Perc_Pos_Neg = pd.read_excel(r'Data_Files\Positive_Negative_Counts_PerGene.xlsx')

Numb_Pos_Int = sum(Perc_Pos_Neg['Number_of_Positive_Interactions_DiauxicShift']) + \
                sum(Perc_Pos_Neg['Number_of_Positive_Interactions_Polarity'])

Numb_Neg_Int = sum(Perc_Pos_Neg['Number_of_Negative_Interactions_DiauxicShift']) + \
                sum(Perc_Pos_Neg['Number_of_Negative_Interactions_Polarity'])

Total_Number_Interactions = Numb_Pos_Int+Numb_Neg_Int

Perc_Pos_Int = Numb_Pos_Int/Total_Number_Interactions
Perc_Neg_Int = Numb_Neg_Int/Total_Number_Interactions

labels = 'Negative', 'Positive'
sizes = [Perc_Neg_Int,Perc_Pos_Int]
cs2 = 'tomato','limegreen'
explode = (0, 0)
fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', textprops={'fontsize': 26},
        shadow=False, startangle=90, colors = cs2 )
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.tight_layout()
