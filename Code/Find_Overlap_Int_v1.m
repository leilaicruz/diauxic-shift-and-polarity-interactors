%% Extract genes that have a genetic interaction with genes related to the Diauxic Shift (DS)
% First, the Diauxic Shift genes and their interactors are imported into a
% table. The interactors are then selected for the type of interaction they
% have with the gene (genetic and NOT physical). The variable
% DS_Interac_Genetic contains all information on the genetic interactors of
% the Diauxic Shift genes. The variable DS_Interac_Genes contains all
% unique genes that have a genetic interaction with Diauxic Shift genes.

% Read in the Diauxic Shift genes and their interactors as a table
DS_Interac = readtable('Genes_Diauxic_Shift_Description_Interactors.xlsx');

% Extract the interaction type (genetic or physical) and select only the
% genetic interactors
DS_Int_Type = DS_Interac.Gene_Interactions_Details_RelationshipType;
DS_Interac_Genetic = DS_Interac(cellfun(@isequal,DS_Int_Type,repmat({'genetic'},size(DS_Int_Type,1),1)),:);

% Get the interactors and store them in a separate variable (for use later
% on in the code)
DS_Int_P2_SecI = DS_Interac_Genetic.Gene_Interactions_Participant2_SecondaryIdentifier;

% Store a table containing the names of the unique genes that have a
% genetic interaction with the diaxic shift genes
DS_Interac_Genes = DS_Interac_Genetic(:,'Gene_Interactions_Participant2_SecondaryIdentifier');
DS_Interac_Genes = unique(DS_Interac_Genes);

%% Extract genes that have a genetic interaction with genes related to Polarity Establishment (Pol)
% First, the Polarity genes and their interactors are imported into a
% table. The interactors are then selected for the type of interaction they
% have with the gene (genetic and NOT physical). The variable
% Pol_Interac_Genetic contains all information on the genetic interactors of
% the Polarity genes. The variable Pol_Interac_Genes contains all
% unique genes that have a genetic interaction with Polarity genes.

% Read in the Polarity genes and their interactors as a table
Pol_Interac = readtable('Polarity_Genes_Description_Interactors.xlsx');

% Extract the interaction type (genetic or physical) and select only the
% genetic interactors
Pol_Int_Type = Pol_Interac.Gene_Interactions_Details_RelationshipType;
Pol_Interac_Genetic = Pol_Interac(cellfun(@isequal,Pol_Int_Type,repmat({'genetic'},size(Pol_Int_Type,1),1)),:);

% Get the interactors and store them in a separate variable (for use later
% on in the code)
Pol_Int_P2_SecI = Pol_Interac_Genetic.Gene_Interactions_Participant2_SecondaryIdentifier;

% Store a table containing the names of the unique genes that have a
% genetic interaction with the diaxic shift genes
Pol_Interac_Genes = Pol_Interac_Genetic(:,'Gene_Interactions_Participant2_SecondaryIdentifier');
Pol_Interac_Genes = unique(Pol_Interac_Genes);

%% Determine the type of interaction of the genes that interact with both 
% This is to find the interactors of both the Diauxic Shift and the
% Polarity genes that overlap, such that we find the genes that have a
% genetic interaction with both (and have possible pleitropic effects with
% regard to these traits)

Overlap_Int = table2array(intersect(Pol_Interac_Genes,DS_Interac_Genes));

%% Generate the matrix that describes the interacting genes
% This part of code loops over the different genes that have a genetic
% interaction with BOTH the Diauxic Shift module and the Polarity module.
% This matrix is built from the perspective of the interactor, which means
% that the column name 'Gene' refers to a gene that interacts with both the
% Diauxic Shift and the Polarity module. As a consequence, each
% 'Interactor' can be assigned to the Diauxic Shift or Polarity module or
% both. Each 'Gene' will be checked if it can be assigned to the Diauxic
% Shift or Polarity module, but might belong to another module (referred to
% as 'Unknown').

% Define the number of fields that contain information about the 'Gene' and
% the 'Interactor' and preallocate the matrix that holds this information
Num_of_Fields = 5;
Gene_Interaction_List = table('Size',[0 Num_of_Fields],'VariableNames',{'Gene','Interactor','Interaction_Type','Interactor_Module','Gene_Module'},'VariableTypes',["string","string","string","string","string"]);

% Loop over all the genes and determine which categories apply to them
for j = 1:size(Overlap_Int)
    
    % Extract the query gene from the list of genes that interact with both
    % modules and determine which Diuxic Shift and Polarity interactors it
    % has (where these interectors occur in the original matrix)
    Query_Gene = Overlap_Int{j,1};
    Query_DS = cellfun(@isequal,DS_Int_P2_SecI,repmat({Query_Gene},size(DS_Int_P2_SecI,1),1));
    Query_Pol = cellfun(@isequal,Pol_Int_P2_SecI,repmat({Query_Gene},size(Pol_Int_P2_SecI,1),1));
    
    % Determine if the the query gene itself belongs to either the Diauxic
    % Shift or the Polarity module
    Query_in_DS_Module = any(cellfun(@isequal,DS_Interac.Gene_SystematicName,repmat({Query_Gene},size(DS_Interac,1),1)));
    Query_in_Pol_Module = any(cellfun(@isequal,Pol_Interac.Gene_SystematicName,repmat({Query_Gene},size(Pol_Interac,1),1)));
    
    % Specify the variables that will be extracted from the original matrix
    % about the information of the query gene and its interactors 
    vars = {'Gene_Interactions_Participant2_SecondaryIdentifier','Gene_SystematicName','Gene_Interactions_Details_Experiment_InteractionDetectionMethod'};
    
    % Extract the information from the original matrix and specify that the
    % interactor is from the Diauxic Shift module
    DS_Info_Query = DS_Interac_Genetic(Query_DS,vars);
    DS_Info_Query.Interactor_Module = repmat({'Diauxic Shift'},size(DS_Info_Query,1),1);
    
    % Do the same for the interactors from the Polarity module
    Pol_Info_Query = Pol_Interac_Genetic(Query_Pol,vars);
    Pol_Info_Query.Interactor_Module = repmat({'Polarity'},size(Pol_Info_Query,1),1);
    
    % Specify if the query gene itself happens to be from the Diauxic Shift
    % or Polarity module, or both.
    if  Query_in_DS_Module && Query_in_Pol_Module 
        
        DS_Info_Query.Gene_Module = repmat({'Diauxic Shift and Polarity'},size(DS_Info_Query,1),1);
        Pol_Info_Query.Gene_Module = repmat({'Diauxic Shift and Polarity'},size(Pol_Info_Query,1),1);
        
    elseif Query_in_DS_Module && ~Query_in_Pol_Module
        
        DS_Info_Query.Gene_Module = repmat({'Diauxic Shift'},size(DS_Info_Query,1),1);
        Pol_Info_Query.Gene_Module = repmat({'Diauxic Shift'},size(Pol_Info_Query,1),1);
        
    elseif ~Query_in_DS_Module && Query_in_Pol_Module
        
        DS_Info_Query.Gene_Module = repmat({'Polarity'},size(DS_Info_Query,1),1);
        Pol_Info_Query.Gene_Module = repmat({'Polarity'},size(Pol_Info_Query,1),1);
        
        
    elseif ~Query_in_DS_Module && ~Query_in_Pol_Module
        
        DS_Info_Query.Gene_Module = repmat({'Unknown'},size(DS_Info_Query,1),1);
        Pol_Info_Query.Gene_Module = repmat({'Unknown'},size(Pol_Info_Query,1),1);
        
    end
    
    % Store the information in the table 'Gene_Interaction_List' that is
    % built as we loop over the different query genes
    Gene_Interaction_List(end+1:end+size(DS_Info_Query,1),1:Num_of_Fields) = DS_Info_Query;
    Gene_Interaction_List(end+1:end+size(Pol_Info_Query,1),1:Num_of_Fields) = Pol_Info_Query;

end

% %% Obtain the different types of genetic interactions that occur in the final matrix
% Interac_Ident = unique(Gene_Interaction_List(:,"Interaction_Type"));
% 
% %% Negative genetic interactions
% Dosage_Growth_Def = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{1,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Dosage_Growth_Def) = "Negative Genetic";
% 
% Dosage_Lethal = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{2,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Dosage_Lethal) = "Negative Genetic";
% 
% Phen_Enhanc = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{5,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Phen_Enhanc) = "Negative Genetic";
% 
% Syn_Growth_Def = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{8,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Syn_Growth_Def) = "Negative Genetic";
% 
% Syn_HaploInsuff = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{9,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Syn_HaploInsuff) = "Negative Genetic";
% 
% Syn_Lethal = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{10,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Syn_Lethal) = "Negative Genetic";
% 
% %% Positive genetic interactions
% Dosage_Rescue = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{3,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Dosage_Rescue) = "Positive Genetic";
% 
% Phen_Supp = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{6,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Phen_Supp) = "Positive Genetic";
% 
% Syn_Rescue = cellfun(@isequal,Gene_Interaction_List.Interaction_Type,repmat(Interac_Ident{11,"Interaction_Type"}, size(Gene_Interaction_List,1),1));
% Gene_Interaction_List.Interaction_Type(Syn_Rescue) = "Positive Genetic";

%% Save the table containing the positive/negative genetic interactions in an excel sheet


