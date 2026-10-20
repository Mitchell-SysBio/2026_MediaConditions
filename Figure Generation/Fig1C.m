%% Figure 1C
% SFS 09-12-25

%% Load variables from excel sheet
% Load file with MGC medians and drug information
fileName =  'FolicAcidDrugScreenMGCs.xlsx';

% Extract drug names
drugNamesTemp = table2array(readtable(fileName, 'Sheet','Plate map key', 'Range', 'B2:B9'));

% Extract condition names and MGC medians
conditionNamesTemp = readtable(fileName, 'Sheet','Median MGCs').Properties.VariableNames(2:12);
MGC_MediansTemp = table2array(readtable(fileName, 'Sheet','Median MGCs', 'Range','C3:M9'));

% Clear variables not necessary anymore
clear fileName

%% Transform to match desired drug class order
% Desired drug order - sorted by class then alphabetically
desiredDrugOrder = {'SD', 'SMM', 'SMX', 'SPY', 'TMP', 'TMP/SMX', 'PYR'};

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

% Remove PYR (was not toxic)
drugNames = drugNames(1:6);
MGC_Medians_ReOrdered = MGC_Medians_ReOrdered([1:6],:);

% Clear variables not necessary anymore
clear count; clear j; clear i; clear drugCatIdx; clear drugNamesTemp; clear MGC_MediansTemp; clear tempIdx; clear desiredDrugOrder; clear drugIdx; clear drugIdx2;
%% Transform again by media order to match SMX increasing resistance order
% Same order as in Fig 1 so importing it
load ../'BMD Dataset'/MGCFig1Variables.mat conditionNames
conditionNamesFig1 = conditionNames;
clear conditionNames

% Reorder to match increasing SMX order from fig 1
% Pull out subset of SMX order from fig 1
[conditionNamesSMXOrder, SMX_order_idx] = intersect(conditionNamesFig1, conditionNamesTemp, 'stable'); % will use idx to pull out colors from mediaColors

% Reordering index
[~, SMX_Idx] = ismember(conditionNamesSMXOrder, conditionNamesTemp);

% Re-sort MGC median by SMX order
MGC_Medians = MGC_Medians_ReOrdered(:,SMX_Idx);
conditionNames = conditionNamesTemp(SMX_Idx);

% Clear variables not necessary anymore
clear MGC_Medians_ReOrdered; clear SMX_Idx; clear conditionNamesTemp; clear conditionNamesFig1; clear conditionNamesSMXOrder;
%% Visualization variables
% Obtaining color/shape subsets from figure 1
load MGCFig1Variables.mat mediaTypesIdx mediaColors shapeOpt
FAScreen_colors = mediaColors(SMX_order_idx,:);
FAScreen_shapeIdx = mediaTypesIdx(SMX_order_idx);

clear mediaColors; clear mediaTypesIdx; clear SMX_order_idx
%% Figure Creation
f1c = figure; hold on;
for i = 1:length(conditionNames)
    x = 1:length(drugNames); % all drugs on x axis
    y = MGC_Medians(:,i);

    % Loop through y values to find if there are two medias with the same MGC for the drug
    for j = 1:length(y)
        if sum(round(y(j),2)==round(MGC_Medians(j,:),2))>1 && sum(round(y(j),2)==round(MGC_Medians(j,:),2)) <3 % if they're equal to any other media for the drug (>1 because it counts itself)
            x(j) = x(j)+(-0.1 + (rand(1)*0.25)); % adding random space that can be positive or negative
        elseif sum(round(y(j),2)==round(MGC_Medians(j,:),2))>=3 %if there's more than two - so they can be more spaced out
            x(j) = x(j)+(-0.25 + (rand(1)*0.45));
        else % if not
        end

    % Then use loop to loop through y values to give them different markers
        p = plot(x(j),y(j), 'LineWidth', 1.5, 'MarkerSize', 10, 'Marker', shapeOpt{FAScreen_shapeIdx(i)}); hold on;
        set(p, 'Color', FAScreen_colors(i,:));
        p.LineStyle = "none";
    end
end

% Miscellaneous figure parameters
xlabel("drugs");
xticks(1:length(drugNames));
xticklabels(drugNames);
ylabel("MGC median");
grid on;
yscale log;
set(f1c, 'Position', [0 0 800 500]);

% Saving in a number of ways
savefig(f1c, ['f1c ' char(date) '.fig']);
%set(f1c, 'Renderer', 'painters');
%saveas(f1c, ['f1c ' char(date) '.svg']);
%exportgraphics(f1c, 'f1c.eps', 'ContentType', 'vector');