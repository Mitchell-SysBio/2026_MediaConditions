%% Supplemental Figure 5B - MGC comparison between WT and eut KOs

%% Load variables
% Load BMD values
MGCsTable = readtable("BMD_eutKOs.xlsx","Sheet","MGCs");
MGCsCell = table2cell(MGCsTable(:,2:end));
MGCsMat = cell2mat(MGCsCell);
clear MGCsTable MGCsCell

% Other variables
koStrains = {'eutI', 'eutN', 'eutQ'};
mediaStarting = {'MAG', 'DMEM', 'MAGV', 'MAGVN', 'MAGN', 'MAGNa'};
mediaDesOrder = {'MAG', 'MAGV', 'MAGN', 'DMEM', 'MAGVN', 'MAGNa'};

%% Prepare data
% Calculate medians and std dev
medianMGCs = []; stdevMGCs = [];
for iStrain = 1:4
    for iCond = 1:length(mediaDesOrder)
        medianMGCs(iStrain,iCond) = median(MGCsMat(iStrain,(1+(4*(iCond-1))):(4+(4*(iCond-1)))));
        stdevMGCs(iStrain,iCond) = std(MGCsMat(iStrain,(1+(4*(iCond-1))):(4+(4*(iCond-1)))));
    end
end

% Reorder to match desired media order
[~,mediaIdx] = ismember(mediaDesOrder, mediaStarting);
newMediaOrder = mediaStarting(mediaIdx);
medianMGCs_Reordered = medianMGCs(:,mediaIdx);
stdMGCs_Reordered = stdevMGCs(:,mediaIdx);

% Divide each KO by WT to get ratio 
mediansDiv = []; stdevDiv = [];
for iKO = 1:length(koStrains)
    mediansDiv(iKO,:) = medianMGCs_Reordered(1+iKO,:)./medianMGCs_Reordered(1,:); % divide by WT
    stdevDiv(iKO,:) = stdMGCs_Reordered(1+iKO,:)./medianMGCs_Reordered(1,:); % divide stdev by the same value as the MGCs because of rule of linear transformations
end

%% Bar plots
fBar = figure; hold on; % creating figure
for iKO = 1:3 % KOs being used
    subplot(1,3,iKO)
    b = bar(mediansDiv(iKO,1:4)); hold on; % create bar plot of medians
    %b.BaseValue = 1;
    errorbar(1:4,mediansDiv(iKO,1:4),stdevDiv(iKO,1:4)); hold on;
    yline(1, '--');
    ylim([0,2.5]);
    xlim([0 5]);
    xticklabels(newMediaOrder(1:4));
    title(koStrains(iKO));
end