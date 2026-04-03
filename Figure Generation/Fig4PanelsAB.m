%% Plot Fig 4 Panels A and B
% SFS 01-25-26
%% load data
load ../'Dose Response Dataset'/doseResponseVar1.mat

%% Panel A - Plot growth curves
% Initial variables
x = time;
numDrugDil = length(avAUC)./length(drugNames); % 12

% Plotting all media conditions by rows, SMX dilutions by columns
f1 = figure; hold on;
tiledlayout(1,12)
plotIdx = smxIdx:numDrugDil; % taking SMX location, and its 12 dilutions
for iDilution = 1:numDrugDil
    nexttile
    p = plot(x, avGC{1}(:,plotIdx(iDilution)));
    grid on; box on;
    ylim([0 0.6]);
end

% Save figure
set(f1, 'position', [295 160 12000 250]);
set(f1, 'Renderer', 'painters');
saveas(f1, 'F4A.fig');
saveas(f1, 'F4A.svg');

%% Calculate IC50s and plot dose response curves - Fig 4 B
% Asssign function variables
tfPlot = 1;
requestedICs = 50;

% Make subfolder
tempDir = 'F4Bplots';
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end 

% Setting up IC50 variable
IC50 = cell(9,1);
for i = 1:length(IC50)
    IC50{i} = zeros(1,2);
end

% Setting up MG control
noDrugMG = avAUC(1, numDrugDil-1);
AUCMG = avAUC(1,1:numDrugDil-1);
normAbsMG = AUCMG./noDrugMG;
MGdrugConc = smxConcMap(1,:);

% Using here with AUC
% First normalizing absorbance to no drug
% Need to loop through each condition for each plate and divide by no drug
for i = 1:length(avGC)
    noDrug = avAUC(i,numDrugDil-1); % subtracting one because 12 is often poor growth
    AUCTemp = avAUC(i,1:numDrugDil-1);
    normAbsorbance = AUCTemp./noDrug;
    [coeffs,aICarray,H2, xpoints, ypoints] = calcDoseResponse(smxConcMap(i,:),normAbsorbance,tfPlot,requestedICs, condsUnique{i}, normAbsMG, MGdrugConc, tempDir); % running calcDoseResponse and getting all of its outputs
    IC50{i}(1,2) = aICarray(2); % IC50 # 2 is abs IC50
end