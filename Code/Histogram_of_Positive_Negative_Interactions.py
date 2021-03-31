# -*- coding: utf-8 -*-
"""
This script calculates the fraction of genetic interactions that are 
positive in the group of shared interactors with at least one positive interaction.

This calculation requires the following data file: 
    
    Genes_With_Positive_Interaction.xlsx
    
"""

import pandas as pd
import scipy 
import seaborn as sns

# Import the datafile (Genes_With_Positive_Interactions) that contains the 
# number of positive and negative interactions per gene for each module 
# (Polarity and Diauxic Shift)
Interac_Types_Count = pd.read_excel(r'Data_Files\Genes_With_Positive_Interaction.xlsx')

# Compare the total number of positive interactions (so regardless of with 
# which module the interactions occur) to the total number of interactions
Interac_Types_Count['Total_Interactions'] = Interac_Types_Count.sum(axis=1)
Interac_Types_Count['Total_Positive_Interact'] = Interac_Types_Count.Positive_DS+Interac_Types_Count.Positive_Pol
Interac_Types_Count['Fraction_Positive'] = Interac_Types_Count.Total_Positive_Interact/Interac_Types_Count.Total_Interactions

# Plot the total number of interactions in a histogram to see if it follows a 
# power-law distribution
sns.set_style('whitegrid')
spacing = 0.05
n_bins = int(1//spacing)

plot = sns.distplot(Interac_Types_Count.Total_Interactions, kde=False, norm_hist=True)

# Plot the fraction of total positive interactions/total interactions as a 
# histrogram and fit it with a lognormal distribution
spacing = 0.05
n_bins = int(1//spacing)

plot = sns.distplot(Interac_Types_Count.Fraction_Positive, kde=False, fit=scipy.stats.lognorm,norm_hist =True, bins = n_bins)

params = scipy.stats.lognorm.fit(Interac_Types_Count.Fraction_Positive)
pdf = lambda x: lognorm.pdf(xvals, *params)

plot.set_xlabel("Fraction positive iteractions",fontsize=14)
plot.set_ylabel("# of genes (scaled)",fontsize=14) 
plot.set(xlim=(0, 1.1))
plot

print(params)