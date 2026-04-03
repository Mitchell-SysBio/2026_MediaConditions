%% Script to produce panels for figure 3
% SFS 01-05-26
% Produces figure panels A, B and C for figure 3
% Figure A: histogram of p values for norm LFC vs potency
% Figure B: visualization of GSEA for p val vs norm LFC results
% Figure C: isolated examples of norm LFC vs potency plots
%% Genes of interest - manual input (for panel C)
genesPanelC = {'metF', 'metH', 'btuB', 'eutQ', 'eutI', 'cchA', 'cchB', 'ypfE', 'panB', 'yadE', 'ydbD'};
geneExampleA = 'ydhR';

%% Load data
% Load genetic screen data
load ../'Genetic Screens Dataset'/GeneticScreensDataset.mat
conditionNames = replace(conditionNames,"_", " ");
clear normLFC_noNaN normLFC LFC padj

% Load correlations between LFC and potency
load ../corrVariables.mat
SMXCond = extractBefore(conditionNames(inxSMX), ' ');
clear conditionNames

% Load pwFinVar.mat
load ../pwFinVars.mat

%% Panel A
% Produce histogram of p values of correlations between KO LFC and WT potency
edgesP = [0:0.1:4];
f1 = figure; hold on;
h2 = histogram(-log10(p), edgesP, 'Normalization','probability');
title('p val');
ylabel('proportion');
xlabel('-log10(p)');
box on; grid on;
set(gcf,'position',[402   600   1100   350]);

% Save figure
saveas(f1, 'fig3a.fig');
% saveas(f1, 'fig3a.png');
% saveas(f1, 'fig3a.svg');

%% Panel B
% Produce circle plot of GSEA results summarized by REVIGO
f2 = figure; hold on;
for iPW =1:height(resP01_10_Sorted)
    curSetSize = pUnivCount_Sorted(iPW); % circle size determined by number of genes present in pw
    plot(iPW, 1,'o','markersize',curSetSize*1.9,'markerfacecolor',[1 1 1],'markeredgecolor','k');
    hold on;
end
grid on; box on;
yticks([0:0.5:2]);
ylim([0 2]);
xticks([1:height(resP01_10_Sorted)]);
xlim([0 height(resP01_10_Sorted)+1]);
xticklabels([resP01_10_Sorted.Name]);
xlabel('pathway');
title('by p val 0.01');
set(gcf,'position',[402   600   1200   290]);

% Save figure
saveas(f2, 'fig3B', 'fig');
saveas(f2, 'fig3B', 'png');
saveas(f2, 'fig3B', 'svg');

%% Panel C
% Producing individual examples of genes of interest LFC vs WT potency

% Generate genesPanelC idx
for iGene = 1:length(genesPanelC) 
    genesCIdx(iGene) = find(contains(strains, genesPanelC(iGene)));
end

% Make subfolder
if ~exist('PanelCPlots/', 'dir') 
    mkdir('PanelCPlots/');
end

% Generate plots for each gene of interest
for iGene = 1:length(genesPanelC)
    fPlot = figure; hold on;

    x = stdPotSMX'; % standardized potency of WT - always the same for each KO
    y = SMX_nLFC(genesCIdx(iGene),:)'; % normalized LFC in SMX for KO

    % Set colors by media
    for iCond = 1:length(SMXCond)
        p1 = plot(x(iCond), y(iCond), '.', 'Color', mediaColors(iCond,:)); hold on;
        p1.MarkerSize = 55;
    end

    % Making fit lines
    % Fit parameters
    lFit = fit(x, y,'poly1'); % making fit line with only non-NaN
    % [rho, p] = corrcoef(x, y);
    rho_fin = r(genesCIdx(iGene));
    p_fin = p(genesCIdx(iGene));

    % Now plotting the fit
    p1 = plot(lFit);
    p1.LineWidth = 0.5;
    p1.Color = 'black';

    str = {['rho = ' num2str(round(rho_fin, 3))],['p = ' num2str(round(p_fin, 3))]}; % labeling with text
    text(-1.5, max(y)-abs(max(y))*0.2, str);

    grid on; box on;
    ylabel('LFC');
    xlabel('potency');
    title(genesPanelC{iGene});

    % putting legend in useful location
    % lgd = legend(SMXCond);
    % if rho_fin >0 % is it positive direction
    %     lgd.Location = "southeast";
    % else % is it negative direction
    %     lgd.Location = "northeast";
    % end

    % Saving plot
    fileName = fullfile('PanelCPlots/', [genesPanelC{iGene}]);
    saveas(fPlot, fileName); % matlab figure
    %saveas(fPlot, fileName, 'svg');
end

%% Generate plot for poor correlation (example for Panel A)

% Generate example gene idex - now just one gene
genesExampleIdx= find(contains(strains, geneExampleA));

% Generate plot
    fPlot = figure; hold on;

    x = stdPotSMX';
    y = SMX_nLFC(genesExampleIdx,:)';

    for iCond = 1:length(SMXCond)
        p1 = plot(x(iCond), y(iCond), '.', 'Color', mediaColors(iCond,:)); hold on;
        p1.MarkerSize = 55;
    end

    % Making fit lines
    % Fit parameters
    lFit = fit(x, y,'poly1'); % making fit line with only non-NaN
    [rho, p] = corrcoef(x, y);
    rho_fin = rho(1,2);
    p_fin = p(1,2);
    %rho_fin = r(genesCIdx(iGene));
    %p_fin = p(genesCIdx(iGene));

    % Now plotting the fit
    p1 = plot(lFit);
    p1.LineWidth = 0.5;
    p1.Color = 'black';

    str = {['rho = ' num2str(round(rho_fin, 3))],['p = ' num2str(round(p_fin, 3))]}; % labeling with text
    text(-1.5, max(y)-abs(max(y))*0.2, str);

    grid on; box on;
    ylabel('LFC');
    xlabel('potency');
    title(geneExampleA);

    % putting legend in useful location
    % lgd = legend(SMXCond);
    % if rho_fin >0 % is it positive direction
    %     lgd.Location = "southeast";
    % else % is it negative direction
    %     lgd.Location = "northeast";
    % end

    % Saving plot
    fileName = fullfile('Example');
    saveas(fPlot, fileName); % matlab figure
    %saveas(fPlot, fileName, 'svg');