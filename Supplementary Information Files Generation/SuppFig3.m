%% Figure S3 - Growth vs Potency
% Load variables
load './SuppFig2ACalculatingGrowthRate/GRVar.mat'
load StdPotFig1Variables.mat
load '../BMD Dataset/MGCFig1Variables.mat' % to get mediaColors

% Clarifying names
stdPot = Z; clear Z;

%% Need to reorder GR order to be in correct order
[~, mediaIdx] = ismember(conditionNames, condsNew);

% Reordering
GRFin = GR(mediaIdx);

% Checking media
% mediaCheck = condsNew(mediaIdx); clear mediaCheck;

%% Remove AU condition
% AU removed to lack of growth that led to poor fit and inaccurate GR
% measurement
AU_idx = find(strcmp('AU', conditionNames));
stdPot(:,AU_idx) = [];
conditionNames(AU_idx) = [];
mediaColors(AU_idx,:) = [];
GRFin(AU_idx) = [];

%% Plotting Fig S1
fGR = figure; hold on;
for iDrug = 1:length(drugNames)
    subplot(3,4,iDrug)
    for iCond = 1:length(GRFin)
        scatter(GRFin(iCond), stdPot(iDrug,iCond), 50, 'LineWidth', 1.5,'MarkerEdgeColor', mediaColors(iCond,:), 'Marker', shapeOpt{mediaTypesIdx(iCond)}); hold on;
    end

    x = GRFin;
    y = stdPot(iDrug,:)'; 
    include =~ isnan(y); % indexing where there aren't NaN values b/c fit can't handle NaN

    % Fit parameters
    lFit = fit(x(include), y(include),'poly1'); % making fit line with only non-NaN
    [rhoEnd, pEnd] = corrcoef(x(include), y(include));
    rho_fin = rhoEnd(1,2);
    p_fin = pEnd(1,2);
    
    % Now plotting the fit
    p = plot(lFit);
    p.LineWidth = 0.5;
    p.Color = 'black';

    str = {['rho = ' num2str(round(rho_fin,3))],['p = ' num2str(round(p_fin,3))]}; % labeling with text
    text(max(x)-0.75, max(y)-0.4, str);

    % Miscellaneous plot parameters
    grid on;
    ylim([floor(min(y))-.5 floor(max(y))+1]);
    yticks(floor(min(y))-.5:0.5:floor(max(y))+1);
    xlim([min(x) max(x)]);
    title(drugNames{iDrug});
    legend('hide');
    xlabel('growth rate (1/hr)');
    ylabel('standardized potency');
    xticks(0:0.25:max(x));
    box on;
end

% Save figure
set(fGR, 'Renderer', 'painters');
savefig(fGR, ['supp3 ' char(date) '.fig']);
%saveas(fGrEndPt, 'supp3', 'svg');