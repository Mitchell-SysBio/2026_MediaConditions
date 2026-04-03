%% Fig 4 Panel C EA experiment
% SFS 01-25-26
%% Load data
load ../'Dose Response Dataset'/doseResponseEAVar1.mat

%% Dose response curves/IC50 calc
% Initial variables
tfPlot = 1;
requestedICs = 50;

% Column count (losing last col due to messiness from evaporation)
colCount = size(drugMap,2)-1; % 11

% Make subfolder
tempDir = 'F4Cplots';
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end 

% Setting up IC50 variable
IC50 = cell(8,1);
for i = 1:length(IC50)
    IC50{i} = zeros(1,2);
end

% Index from experiment with no EA for SMX or TMP
SMXIdx = 4;
TMPIdx = 8;

% Using here with AUC
% First normalizing absorbance to no drug
% Need to loop through each condition for each plate and divide by no drug
for i = 1:size(avAUC,1)
    if i <= 4 % SMX conditions
        noDrugCont = avAUC(SMXIdx, colCount); % 4 for SMX no EA
        AUCCont = avAUC(SMXIdx,1:colCount);
        normAbsCont = AUCCont./noDrugCont;
        contDrugConc = drugMap(SMXIdx,1:colCount);
    else % TMP conditions
        noDrugCont = avAUC(8, colCount); % 8 for TMP
        AUCCont = avAUC(8,1:colCount);
        normAbsCont = AUCCont./noDrugCont;
        contDrugConc = drugMap(8,1:colCount);
    end
    noDrug = avAUC(i,colCount); 
    AUCTemp = avAUC(i,1:colCount);
    normAbsorbance = AUCTemp./noDrug;
    [coeffs,aICarray,H2, xpoints, ypoints] = calcDoseResponse(drugMap(i,1:colCount),normAbsorbance,tfPlot,requestedICs, condMap{i}, normAbsCont, contDrugConc, tempDir); % running calcDoseResponse and getting all of its outputs
    IC50{i}(1,2) = aICarray(2); % IC50 # 2 is abs IC50
end