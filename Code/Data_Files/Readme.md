## Data files origin 

- `Overlap_First_Interactors_NoDuplicates.xlsx` 
    - script that use it: `Diversity_Histogram.py`
    - source:
        - what is described in this file? 

This file contains 4 columns ('Gene', 'Interactor','Interaction_Type' and 'Interactor_Module')
and describes the genes that interact with both modules (polarity+diauxic shift).
These genes are contained in the column 'Gene', and their interaction partner from either of the two modules is contained 
in the column 'Interactor'.

Because in some cases a polarity gene will directly interact with a diauxic shift gene (or vice versa), 
these interactions would appear twice in this dataset (once with the polarity gene as 'Gene' and the diauxic shift gene as 'Interactor' 
and once with their roles reversed). This dataset has been cleared of these to avoid double counting. 
The dataset Overlap_First_Interactors contains the original dataset with duplicates

    - how was created this file? 
    