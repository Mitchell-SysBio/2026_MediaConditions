%% Genetic screens screen initial data preparation script
% SFS 11-19-2025

%% User definitions
cutoffLFC = log2(1.5);
cutoffPadj = 0.25;
%% File locations
dataFolder = 'all_barcodes/';
dataList = dir([dataFolder '*.csv']);
nConditions = length(dataList);

%% Read files and extract corrected counts, LFC and padj
corrCounts = zeros(3680,nConditions); LFC = zeros(3680,nConditions); % set the number of expected strains based on library size
padj = zeros(3680,nConditions); conditionNames = cell(nConditions,1); % additional empty variables
for iFile =1:nConditions % loop through conditions to get LFC
    curName = dataList(iFile).name;
    curT = readtable([dataFolder curName]); % read screen condition with all knockout LFCs
    splitName = split(curName,'_LFC.csv'); % obtain condition name
    if height(curT) ~= 3680 % check for correct number of strains in library
        error(['Error: expecting data for 3680 strains but got ' ...
            num2str(height(curT)) ' for ' splitName{1} '. Check R output'])
    end
    corrCounts(:,iFile) = curT.baseMean;
    LFC(:,iFile) = curT.log2FoldChange; % collect log2fold change for each knockout
    padj(:,iFile) =  curT.padj; % collect p adj value for each log2fold change
    splitName = split(splitName(1),'ASKA_');
    conditionNames(iFile,1) = splitName(2); % save unique condition name
end
strains = curT.Var1;
clear curT curName splitName corrCounts

% Remove strains with 8+ NaN values!
tfRemove = sum(~isnan(LFC),2)<=12; % check whether 12/20 conditions have non NaN values
strains(tfRemove,:)=[];LFC(tfRemove,:)=[];padj(tfRemove,:)=[]; % remove knockouts with 8+ NaN values across conditions

%% Removing BHI TMP
% BHI TMP removed due to poor quality of sequencing results
% BHI TMP index
BHI_TMP_idx = 2; % manually identifies index for BHI-TMP condition

% deleting that index in all major variables
LFC(:,BHI_TMP_idx) = [];
conditionNames(BHI_TMP_idx) = [];
dataList(BHI_TMP_idx) = [];
nConditions = nConditions - 1;
padj(:,BHI_TMP_idx) = [];

%% Normalizing LFC
stdLFC = std(LFC,'omitmissing'); % standard deviation of LFC
normLFC = LFC./stdLFC; % normalize by dividing by std dev
%normLFC_noNaN = knnimpute(normLFC',10)'; 

%% Media colors variable
% media colors in order of data organization
mediaColors = [111 6 2; % BHI
    255 180 160; %DMEM
    10 130 250; % MAG-Na
    0 170 30; %MAFr
    120 210 250; % MAGly
    200 85 150 % MG
    150 70 150; %MAGN
    173 220 11; % MAG
    255 131 30; % MAGV
    200 6 2; % MAGVN
    ]./265;

%% saving variables
save('GeneticScreensDataset.mat', 'conditionNames', 'LFC', 'normLFC', 'normLFC_noNaN', 'strains', 'padj', 'mediaColors');