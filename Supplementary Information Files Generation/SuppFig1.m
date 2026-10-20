%% Figure S1
% Bar plots showing replicates for BMD screens
% plot replicates for supplemental figure 1

% Load BMD well data
wellCalls1 = table2cell(readtable('../BMD Dataset/MGCRecordsFullBMDScreen.xlsx', 'Sheet','Well calls each plate'));
wellCalls2 = wellCalls1(2:12,:); % removing replicate number in top row
drugNames = wellCalls2(:,1); % isolating drug names
wellCalls = wellCalls2(:,2:64); % removing drug names
clear wellCalls1 wellCalls2

% Reorder drugs to match fig1 order
drugOrder = {'Aztreonam', 'Cefdinir', 'Azithromycin', 'Kanamycin', 'Tetracycline HCl', 'Ciprofloxacin', 'Norfloxacin', 'Sulfamethoxazole', 'Trimethoprim', 'Colistin', 'Nitrofurantoin'};
[~,drugIdx] = ismember(drugOrder, drugNames);
newDrugOrder = drugNames(drugIdx);
wellCallsOrdered = wellCalls(drugIdx,:); % reordering well calls to match new drug order matching fig 1
clear drugNames drugIdx

% Establish drug abbreviations for plot
drugAbbrev = {'ATM', 'CFDN', 'AZM', 'KAN', 'TET', 'CIP', 'NOR', 'SMX', 'TMP', 'CST', 'NIT'};

% Load media names in order
sheet1 = readtable('../BMD Dataset/MGCRecordsFullBMDScreen.xlsx', 'Sheet','Median MGCs');
mediaNames1 = sheet1.Properties.VariableNames;
mediaNames = mediaNames1(2:end);
clear sheet1 mediaNames1

%% Plot

iCond = 1; % starting at condition 1 (BHI)
tempWellCallsInt = zeros(12,3); % setting empty variable

fBar = figure; hold on; % creating figure
for i = 1:length(mediaNames) % looping through all 21 conditions (3 replicates per condition)
    subplot(7,3,i); % plotting each media
    tempReplicates = ((iCond-1)*3)+1:((iCond-1)*3)+3; % isolating the 3 replicates for that condition
    for j = 1:length(drugOrder) % looping through all drugs
        tempWellCallsInt(j,:) = cell2mat(wellCallsOrdered(j,tempReplicates));
    end
    bar(tempWellCallsInt, 'hist');  % create bar plot of well calls
        % xlabel('drugs');
        ylabel('well');
        ylim([0, 9]);
        grid on;
        yticks(0:2:9);
        title(mediaNames(i));
        iCond = iCond + 1; % move to next condition for next loop
        xlim([0.25 11.75]);
        xticklabels(drugAbbrev);
end
sgtitle('bar plots of well of max growth for each replicate')
%savefig(fBar2, 'bar plots comparing replicates.fig');