%% Plot Fig 4 Panel D
% SFS 01-25-26
%% Load data
load ../'Dose Response Dataset'/doseResponseEA2Var.mat

%% Calculate IC50s and plot dose response curves - Fig 4 D
% Initial variables
tfPlot = 1;
requestedICs = 50;

% Make subfolder
tempDir = 'F4Dplots';
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end 

% Column count (losing last col due to messiness from evaporation)
colCount = size(drugMap,2)-1; % 11

% Setting up IC50 variable
IC50 = cell(6,1);
for i = 1:length(IC50)
    IC50{i} = zeros(1,2);
end

% Using here with AUC
% First normalizing absorbance to no drug
% Need to loop through each condition for each plate and divide by no drug
% In this experiment setup every other row is no EA control
% So i+1 (even rows) is no EA control row, and i (odd rows) are 10 uM EA
for i = 1:2:size(avAUC,1)
    % No EA control (i+1)
    noDrugCont = avAUC(i+1, colCount);
    AUCCont = avAUC(i+1,1:colCount);
    normAbsCont = AUCCont./noDrugCont;
    contDrugConc = drugMap(1:colCount);
    
    % EA condition
    noDrug = avAUC(i,colCount);
    AUCTemp = avAUC(i,1:colCount);
    normAbsorbance = AUCTemp./noDrug;
    [coeffs,aICarray,H2, xpoints, ypoints] = calcDoseResponse(drugMap(1:colCount),normAbsorbance,tfPlot,requestedICs, condsUnique{i}, normAbsCont, contDrugConc, tempDir); % running calcDoseResponse and getting all of its outputs
    IC50{i}(1,2) = aICarray(2); % IC50 # 2 is abs IC50
end
