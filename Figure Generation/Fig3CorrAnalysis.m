%% Analysis of growth curve v potency for use in Fig 3
% SFS 01-05-26
% Script does the following: 
% - loads potency and norm LFC data and isolates SMX conditions
% - calculates correlation coefficient and p value for those correlations
% - prepares .txt documents to be used for GSEA analysis using chosen p
% value cutoff
%% Manual definitions
% p val cutoff for potency vs norm LFC data
pCutoff = 0.01;

%% Load datasets
load ../'Genetic Screens Dataset'/GeneticScreensDataset.mat
conditionNamesASKA = replace(conditionNames,"_", " ");
clear normLFC_noNaN conditionNames

% Manual rename of conditions (without drug name) to match potency format
conds = {'BHI', 'DMEM', 'MAG_Na', 'MAFr', 'MAGly', 'MG', 'MAGN', 'MAG', 'MAGV', 'MAGVN'};

%% Load potency data
load ../'BMD Dataset'/StdPotFig1Variables.mat
conditionNamesPot = conditionNames;
stdPot = Z;
clear conditionNames Z

%% Isolate SMX from potency and reorder
%Identify SMx conditions and SMX potency
inxSMXPot = find(contains(drugNames, 'Sulfamethoxazole'));
stdPotSMX = stdPot(8,:);

% Reorder so that they match across datasets
for iCond = 1:length(conds)
    reorderIdx(iCond) = find(strcmp(conditionNamesPot, conds{iCond}));
end
stdPotSMX = stdPotSMX(reorderIdx);

% Double check successful indexing
condCheck = conditionNamesPot(reorderIdx);
if sum(strcmp(condCheck, conds)) < length(conds)
    disp('index ordering failure');
    return;
else
end

clear inxSMXPot reorderIdx iCond condCheck

%% Isolating SMX hits only
% Isolating SMX normLFC
inxSMX = find(contains(conditionNamesASKA, 'SMX'));
SMX_nLFC = normLFC(:,inxSMX);

%% Calculating correlation between std potency and norm LFC
% Correlation calcualted for each ko strain
for iGene = 1:length(strains)
    [r(iGene), p(iGene)] = corr(stdPotSMX', SMX_nLFC(iGene,:)', 'Rows','pairwise');
end

%% Print off datasets to funnel through GSEA
% Universe - all ASKA strains
univStrains = strains;
% Open the file for writing
fileID = fopen('Universe.txt', 'w');

% Write each string to the file, followed by a newline character
for i = 1:length(univStrains)
    fprintf(fileID, '%s\n', univStrains{i});
end

% Close the file
fclose(fileID);

% Hist list based on pval cutoff
% Find hit list
pCutIdx = find(p<=pCutoff);
pCutStrains = univStrains(pCutIdx);

% Open the file for writing
fileID = fopen(['pCorr.txt'], 'w');

% Write each string to the file, followed by a newline character
for j = 1:length(pCutStrains)
    fprintf(fileID, '%s\n', pCutStrains{j});
end

% Close the file
fclose(fileID);

%% Save variables for use in fig 3
save('corrVariables.mat', 'r', 'p', 'SMX_nLFC', 'stdPotSMX', 'inxSMX');