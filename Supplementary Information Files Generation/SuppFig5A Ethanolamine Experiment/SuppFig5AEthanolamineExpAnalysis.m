%% Ethanolamine platereader experiment analysis
% Data is for BHI and M9 vit - only M9 vit data is used
% SFS 10-19-2026

%% Loading data and initial normalizations
% Load the data
dataTable1 = readtable("96Well_GrowthCurve_OD600_NUNC_20250607_110738.xlsx");

% Cutoff time for use in analysis
realVal = 65; % taking the first 65 measurements

time = 0:0.33:0.33*(realVal-1); % manually setting time

% Moving to arrays
dataArrays = {};
for i = 1
    dataArrays{i} = table2array(dataTable1(1:realVal,4:end));
end

% Subtract the mean optical density at the first time point from all other measurements
for i = 1:length(dataArrays)
    mean_OD_T1 = mean(dataArrays{i}(1,49:96));
    max_optical_density_ever(i) = max(max(dataArrays{i}(:,:))); %will also use this value later
    dataArrays{i}(:,:) = dataArrays{i}(:,:) - mean_OD_T1;
end

% Conditions
condNames = {'BHI', 'M9 vit'};

clear mean_OD_T1 
clear i 

%% Calculating area under curve (AUC)
AUC = {};
for i = 1:length(dataArrays)
    AUC{i} = zeros(96,1);
end

for i = 1:length(dataArrays)
    for j = 1:96
        y = dataArrays{i}(:, j);
        AUC{i}(j) = trapz(time,y);
    end
end

%% Taking averages of replicates
% Cycling through media conditions and drug conditions to take averages
for i = 1:2:8 % 4 conditions (BHI/M9 vit and SMX/no drug)
    for j = 1:12 % 12 EA conc.
        averagesTemp(:,j) = mean(dataArrays{1}(:,[(j+(i-1)*12), 12+j+((i-1)*12)]),2);
        averagesAUC(i,j) = mean(AUC{1}([(j+(i-1)*12), 12+j+((i-1)*12)]));
    end
    dataAverages(:,((i-1)*12)+1:((i-1)*12)+12) = averagesTemp;
end

avAUCSMXVit = averagesAUC(5,:);

%% Plotting growth curves
% Compound information
EAConc = [40000./(2.^(0:9)), 0];

% Colors
myTurbo = turbo(length(EAConc));

% Drugs
drugNames = {'SMX', 'no drug'};

% Plot
x = time; % setting times
for iPlate = 1:length(dataArrays)
    for iRow = 5:2:8 % cycling through rows of data averages (ONLY need MAGV condition)
        rowFig = figure; hold on;
        y = dataAverages(:,((iRow-1)*12)+1:((iRow-1)*12)+11);
        for iCol = 1:length(myTurbo) % cycling through rows to change colors
            plot(x,y(:,iCol), '-', 'Color', myTurbo(iCol,:), 'LineWidth', 2); hold on;
        end
        legend([string(EAConc(1:11))], "Location","southeast");
        grid on; box on;
        xlim([0 max(time)]);
        ylim([-0.1 0.8]);
        sgtitle([drugNames{mod((floor((iRow-1)./2)+1)-1,2)+1} ' ' condNames{floor((iRow-1)/4)+1}]);
        %savefig(rowFig, [drugNames{mod((floor((iRow-1)./2)+1)-1,2)+1} ' ' condNames{floor((iRow-1)/4)+1}]);
    end
end

%% Inset plot with AUC values
% EA conc
EAConc = [20, 10, 5, 2.5, 1.25, 0.625, 0.3125, 0.156, 0.078, 0];

% EA conc x axis, av AUC y axis (not including 40 mM conc)
fig2 = figure; hold on;
plot(EAConc, avAUCSMXVit(2:11), LineStyle="none", Marker=".", MarkerSize=40);
xscale("log");
xlim([0 20]);
%xticks([0 1 10 20]);
ylim([0.4 4.5]);
yticks(1:4);
grid on; box on;
ylabel('AUC');
xlabel('log([EA]) (mM)')