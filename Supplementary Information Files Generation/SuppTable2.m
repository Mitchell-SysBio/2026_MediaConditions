%% Supplementary Table 2 Generation
% SFS 01-26-26
%% Load GSEA variables
% using data from Figure 2 GSEA
load gseaResults.mat

% All ontologies
ontologies = gseaOutput.OntologyName;
%% Loop through ontologies to make excel file

% Loop through each gene ontology database
for iOnt = 1:length(ontologies)
    % Setting next ontology database
    databaseDes = ontologies(iOnt);
    ontIdx = find(strcmp(databaseDes, gseaOutput.OntologyName));

    % Establishing gene universe and hits in database
    genesUniverse = gseaOutput.genesUniverse(ontIdx);
    genesUniverse = genesUniverse{1,1};
    genesPresent = gseaOutput.genesPresent(ontIdx);
    genesPresent = genesPresent{1,1};
    pwNames = gseaOutput.PwName(ontIdx);
    pwNames = pwNames{1,1};
    pFDR = gseaOutput.pvalFDR(ontIdx);
    pFDR = pFDR{1,1};

    % Calculating ratio of hit genes in pw / genes in pw universe and creating table
    % Prepare empty table variable to be filled in loop
    gseaT = table(...
        strings(0,1), ...  % PW_name
        zeros(0,1),   ...  % pval_FDR
        strings(0,1), ...  % Gene_Hits
        strings(0,1), ...  % Gene_Universe
        zeros(0,1),   ...  % ratio
        'VariableNames', {'PW_name', 'pval_FDR', 'Gene_Hits', 'Gene_Universe', 'Ratio_Hits_Universe'});

    ratios = [];
    for iPW = 1:length(pwNames) % for each pathway
        tempGenesUniv = genesUniverse{iPW};
        tempGenesPres = genesPresent{iPW};
        overlapCount = sum(ismember(tempGenesUniv, tempGenesPres));
        ratios(iPW) = overlapCount/length(tempGenesUniv);

        % Filling table
        gseaT.PW_name(iPW) = pwNames{iPW};
        gseaT.pval_FDR(iPW) = pFDR(iPW);
        gseaT.Gene_Hits(iPW) = strjoin(tempGenesPres, ' '); % combining to put into table
        gseaT.Gene_Universe(iPW) = strjoin(tempGenesUniv, ' ');
        gseaT.Ratio_Hits_Universe(iPW) = ratios(iPW);
    end

    % Save table to excel file
    filename = 'SuppTables2.xlsx';
    newSheetName = databaseDes;
    writetable(gseaT, filename, 'Sheet',newSheetName{1});

    clearvars -except ontologies iOnt gseaOutput
end