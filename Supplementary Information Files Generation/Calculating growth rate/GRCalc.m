%% Calculating growth rate from tecan experiments
% SFS 12-07-25

% Calculates growth rate from Tecan growth curve experiments 
% using the function fitKishonyRateFunc.m
%% Load variables
% Load the Tecan growth data
dataTable1 = readtable('20240323_results.xlsx');
dataTable2 = readtable('20240404_results.xlsx');
dataTable3 = readtable('20240418_results.xlsx');
dataTable4 = readtable('20240423_results.xlsx');

% All data tables
dataTables = {dataTable1, dataTable2, dataTable3, dataTable4};

% Converting to array and removing excess variables
dataArrays = {};
for i = 1:length(dataTables)
    dataArrays{i} = table2array(dataTables{i}(:,4:end));
end

% Subtract the mean optical density at the first time point from all other measurements
for i = 1:length(dataArrays)
    abs_t1 = mean(dataArrays{i}(2:4,:)); % average of first 3 for noise
    dataArrays{i}(:,:) = dataArrays{i}(:,:) - abs_t1;
    dataArrays{i}(1,:) = 0;
end

% Importing metadata and conditions
conds = table2cell(readtable("gr_conditions.xlsx")); % list of all media conditions

% Metadata sheets with plate maps
for iPlate = 1:length(dataArrays) % for identifying how well was used based on plate maps
    metadata{iPlate} = table2cell(readtable("gr_metadata.xlsx", "Sheet", ['plate ' num2str(iPlate)]));
end

% Final time of experiment run = 20 hours (41 measurements)
finTP = 41;

% Array of times for measuring GR - MANUALLY identified time that was 1 hour after plateued growth
% in hours
finTPGRArray = [7, 18, 10, 8, 10, 19, 9, 15, 11, 12, 11, 11, 13, 11, 11, 11, 8, 9, 9, 11, 8];

%% Convert values from abs to OD
load ODtoCFU120.mat % load calibration curve for 120 uL
for iPlate = 1:length(dataArrays)
    dataArraysOD{iPlate} = (dataArrays{iPlate} - lm.Coefficients.Estimate(1))./lm.Coefficients.Estimate(2);
end

%% Change metadata to be vectorized
% Translate metadata orders to well orders
metadata_vector = cell(1,4);
for iPlate = 1:length(dataArraysOD)
    metadata_vector{iPlate} = cell(96,1);
    currPlate = metadata{iPlate};
    for iRow = 1:8
        metadata_vector{iPlate}(((iRow-1)*12)+1:((iRow-1)*12)+12) = currPlate(iRow,:);
    end
end

%% Calculate mean and std dev (and remove outliers)
% Calculate unique conditions per plate
for iPlate = 1:length(dataArraysOD)
    uniqueConds{iPlate} = unique(metadata{iPlate});
    uniqueConds{iPlate}(strcmp(uniqueConds{iPlate}, 'empty')) = []; % remove empty
end

% Calculate mean and stdev per condition - and removing outliers
for iPlate = 1:length(dataArraysOD)
    for iCond = 1:length(uniqueConds{iPlate})
        tempIdx = find(strcmp(uniqueConds{iPlate}(iCond), metadata_vector{iPlate})); % finding where cond is
        tempWells = dataArraysOD{iPlate}(:,tempIdx); 
        temp_mean = mean(tempWells,2); % temporary mean could still have outliers
        temp_std = std(tempWells,[],2); % temp std could still have outliers
        if temp_std(25,:) >= 0.05 % if std is high at midpoint of growth then there's probably an outlier
            outlier_idx = find(tempWells(25,:) >= temp_std(25,:)+temp_mean(25,:)); % find where it's over 1 std higher
            tempWells(:,outlier_idx) = NaN; % remove outlier by making it NaN
        else % if std isn't high, change nothing
        end
        gc_av{iPlate}(:,iCond) = mean(tempWells(1:finTP,:),2, 'omitNan'); % saving as gc_av with either changed temp or unchanged temp
        gc_std{iPlate}(:,iCond) = std(tempWells(1:finTP,:),[],2, 'omitNan');
    end
end

% Condense to one array
meansArray = []; stdArray = []; condsNew = {};
start = 1;
for iPlate = 1:length(dataArraysOD)
    arrayNum = length(uniqueConds{iPlate});
    meansArray(:,start:start+arrayNum-1) = gc_av{iPlate};
    stdArray(:,start:start+arrayNum-1) = gc_std{iPlate};
    condsNew(start:start+arrayNum-1) = uniqueConds{iPlate};
    start = start+arrayNum;
end

%% Remove duplicate or unwanted conditions
% Some conditions were used multiple times across plates as a control for
% consistency across experiment days. Since they were consistent, their
% replicates can be removed.

% Additionally, condition M9 + gluconate was removed from all further
% experiments due to poor growth

% Identifying conditions that occur multiple times
uniqueAll = unique(condsNew);
for iCond = 1:length(uniqueAll)
    tempCond = strcmp(uniqueAll{iCond}, condsNew);
    condCount(iCond) = sum(tempCond);
end
duplicates = uniqueAll(condCount >= 2);

% Finding duplicate conditions
for iDup = 1:length(duplicates)
    idx = find(strcmp(duplicates{iDup}, condsNew));
    % checking for IGNORE condition - want to remove all cases of IGNORE
    if strcmp(duplicates{iDup}, 'IGNORE')
        removeIdx = idx;
    else % only want to remove duplicates
        removeIdx = idx(2:end);
    end
    % removing wells
    meansArray(:,removeIdx) = [];
    stdArray(:, removeIdx) = [];
    condsNew(removeIdx) = [];
end

clear conds gc_av gc_std 

%% Use fit function to calculate growth rate
% Calculating growth rate
GR = zeros(size(meansArray,2),1);

% Looping through to get get parameters from fit function
% Using new index for ordering
for iCond = 1:size(meansArray,2)
   
    % Isolating condition up to desired time point
    y = meansArray(1:(finTPGRArray(iCond)*2)+1, iCond)'; % *2 because of being in 30 min intervals, +1 because of 0
    x = 0:0.5:finTPGRArray(iCond);
    
    % Using curve fitting function to isolate growthRate
    [sol,~,growthRate, ~, ~, ~] = fitKishonyRateFunc(x, y);
    GR(iCond) = growthRate;
    fits{iCond} = sol;

    % Also isolating AUC of growth curve as alernative metric (optional)
    AUC(iCond) = trapz(y);
end

%% Plotting fit curves (for confirming adequate fit - optional)
% for i = 1:length(fits)
%     subplot(4,6,i)
%     plot(0:0.5:(finTP/2)-0.5, meansArray(1:(finTP),i), 'Marker','.', 'LineStyle','none', 'Color', 'r'); hold on;
%     plot(0:0.5:finTPGRArray(i), meansArray(1:(finTPGRArray(i)*2)+1,i), 'Marker','.', 'LineStyle','none', 'Color', 'k'); hold on;
%     p = plot(fits{i});
%     p.Color = 'k';
%     title(condsNew{i})
%     text(10, 0.25, num2str(GR(i)));
%     xlabel('time (hours)');
%     ylabel('OD at 600 nm');
%     grid on;
% end

%% Save growth rate variable
save('GRVar.mat', 'GR', 'condsNew', 'AUC');