%% Figure S2B
% Comparing endpoint and growth rate 

% Load both datasets and rename for clarity
% Non-shaking growth rate measurement
load StatGrowth.mat
media_stat = media_order;
clear media_order;
% Shaking culture growth rate measurement
load ../SuppFig2ACalculatingGrowthRate/GRVar.mat
media_GR = condsNew;
clear condsNew;
clear AUC;

%% Reorder media orders to match
% Manually rename media for non-shaking media for clarity
media_stat_2 = {'MAG', 'MAGVN', 'MA', 'MG', 'MAFr', 'MAGly', 'MAP', 'MAR', 'MAFu', 'MAGV', 'MAGN', 'MAG_Na', 'MAG_Tr', 'RPMI', 'DMEM', 'AU', 'MHB', 'MHB_S', 'MHB_F', 'BHI', 'LB'};

% Get index for altering media order
[~, idx] = ismember(media_GR, media_stat_2);

% Reorder media according to GR order
media_stat_reordered = media_stat_2(idx);
endPt_reordered = endPt_norm(idx);

%% Plot GR vs endpoint growth
% Colors to distinguish media
myTurbo = turbo(21);

% Plotting
F = figure; hold on;
for i = 1:length(GR)-1 % not AU at the end
    plot(GR(i),endPt_reordered(i), '.', 'MarkerSize', 20, 'Color', myTurbo(i,:)); hold on;
end

    x = GR(1:20);
    y = endPt_reordered(1:20)'; 
    include =~ isnan(y); % indexing where there aren't NaN values b/c fit can't handle NaN

    % Fit parameters
    lFit = fit(x(include), y(include),'poly1'); % making fit line with only non-NaN
    [rhoEnd, pEnd] = corrcoef(x(include), y(include));
    rho_fin = rhoEnd(1,2);
    p_fin = pEnd(1,2);
    
    % Plotting the fit
    p = plot(lFit);
    p.LineWidth = 0.5;
    p.Color = 'black';

    str = {['rho = ' num2str(round(rho_fin,3))],['p = ' num2str(round(p_fin,3))]}; % labeling with text
    text(max(x)-0.75, max(y)-0.4, str);

    grid on; box on;
    xlabel('growth rate shaken culture');
    ylabel('stationary growth absorbance measurement');
