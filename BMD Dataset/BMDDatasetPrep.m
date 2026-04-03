%% BMD Dataset preparation
% Written by SFS 12-09-25
% Initial basic analysis of BMD screen results to be utilized in later
% figure development
%% Load variables from excel sheet
% File with MGC medians and drug information
fileName =  'MGCRecordsFullBMDScreen.xlsx';

% Extract drug names
drugNamesTemp = readtable(fileName, 'Sheet','Plate map key').Properties.VariableNames(2:12);

% Extract condition names and MGC medians
conditionNamesTemp = readtable(fileName, 'Sheet','Median MGCs').Properties.VariableNames(2:22);
MGC_MediansTemp = table2array(readtable(fileName, 'Sheet','Median MGCs', 'Range','C3:W13'));

% Clear variables not necessary anymore
clear fileName
%% Transform to match desired drug class order
% Desired drug order - sorted by class then alphabetically
desiredDrugOrder = {'Aztreonam', 'Cefdinir', 'Azithromycin', 'Kanamycin', 'TetracyclineHCl', 'Ciprofloxacin', 'Norfloxacin', 'Sulfamethoxazole', 'Trimethoprim', 'Colistin', 'Nitrofurantoin'};

% Changing MGC to match this order
[~, drugIdx] = ismember(drugNamesTemp, desiredDrugOrder);

% Using index above to reorder the drug name and the MGC median variables
drugNames = {}; MGC_Medians_ReOrdered = [];

count = 1;
for i = 1:length(drugIdx)
    tempIdx = find(drugIdx == i);
    for j = tempIdx
        drugNames{count} = drugNamesTemp{j}; % drug name reordering
        MGC_Medians_ReOrdered(count,:) = MGC_MediansTemp(j,:); % median MGC reordering
        drugIdx2(count) = drugIdx(j); % drug index updated to new order as sanity check
        count = count + 1;
    end
end

% Clear variables not necessary anymore
clear count; clear j; clear i; clear drugCatIdx; clear drugNamesTemp; clear MGC_MediansTemp; clear tempIdx; clear desiredDrugOrder; clear drugIdx; clear drugIdx2;

%% Transform again by media order to match SMX increasing resistance order
% Isolate SMX conc
MGC_SMX = MGC_Medians_ReOrdered(8,:);

% Identify max order index for SMX
[~, MGC_SMX_Idx] = sort(MGC_SMX, 'descend');

% Re-sort MGC median by SMX order
MGC_Medians = MGC_Medians_ReOrdered(:,MGC_SMX_Idx);
conditionNames = conditionNamesTemp(MGC_SMX_Idx);

% Clear variables not necessary anymore
clear MGC_Medians_ReOrdered; clear MGC_SMX; clear MGC_SMX_Idx; clear conditionNamesTemp;
%% Visual variables
mediaColors = [111 6 2; % reds 5 total
    200 6 2;
    255 20 20;
    255 120 120;
    255 180 160; ...
    255 155 100; % oranges-yellows 6 total
    255 131 30;
    250 180 30;
    225 180 11;
    250 230 55; % 10
    205 230 11;% greens 3 total
    173 220 11;
    0 170 30;
    120 230 200; % blues 8 total
    120 210 250;
    10 130 250;
    30 86 225;
    0 0 150;
    120 90 165; % 3 purple
    150 70 150;
    200 85 150
    ]./265;

% Indexing for changing marker type by media type
mediaTypes = {'defined', 'undefined'};
mediaTypesIdx = [2 1 1 2 1 1 1 2 2 1 2 1 1 1 1 1 1 1 1 1 1];
shapeOpt = {"o", "diamond"}; % defined = circle, undefined = diamond

%% Calculating standardized potency
% Standardize the drug conc data (mean=0, STD =1)
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x, 'omitnan')), std(x, 'omitnan'));
Z = -zscor_xnan(log10(MGC_Medians)')'; % Z is normalized potency

%% Save variables
save("MGCFig1Variables", 'MGC_Medians', 'conditionNames', 'mediaColors', 'mediaTypesIdx', "shapeOpt");
save("StdPotFig1Variables", 'Z', 'conditionNames', 'drugNames');
