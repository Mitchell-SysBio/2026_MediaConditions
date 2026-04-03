%% Supp Fig 2
% SFS 01-252-26
%% Load data
load ../'Dose Response Dataset'/doseResponseSupp2.mat

%% Calculate IC50s for additional drugs
% Initial varaibles
x = time;
numDrugDil = length(avAUC)./length(drugNames); % 12

% Misc variables
tfPlot = 1;
requestedICs = 50;

% New folder
tempDir = 'Supp2plots';
if ~exist(tempDir, 'dir')
    mkdir(tempDir);
end 

% Setting up IC50 variable
IC50 = cell(27,1);
for i = 1:length(IC50)
    IC50{i} = zeros(1,2);
end

% Using AUC to calc IC50
% First normalizing absorbance to no drug
% Need to loop through each condition for each plate and divide by no drug

% Setting up MG control
count = 1;
% Looping through each of the 3 drugs
for iDrug = 1:3
    noDrugMG = avAUC(1, ((iDrug+1)*numDrugDil)-1);
    AUCMG = avAUC(1,(iDrug*numDrugDil)+1:((iDrug+1)*numDrugDil)-1);
    normAbsMG = AUCMG./noDrugMG;
    MGdrugConc = addlConcMaps(iDrug,:);

    % Now setting up conditions with drug
    for i = 1:length(avGC)
        noDrug = avAUC(i,(iDrug+1)*numDrugDil-1); % subtracting one because 12 is often poor growth
        AUCTemp = avAUC(i,(iDrug*numDrugDil)+1:((iDrug+1)*numDrugDil)-1);
        normAbsorbance = AUCTemp./noDrug;
        [coeffs,aICarray,H2, xpoints, ypoints] = calcDoseResponse(addlConcMaps(iDrug,:),normAbsorbance,tfPlot,requestedICs, [condsUnique{i} drugNames{iDrug+1}], normAbsMG, MGdrugConc, tempDir); % running calcDoseResponse and getting all of its outputs
        IC50{count}(1,2) = aICarray(2); % IC50 # 2 is abs IC50
    end
end