%% Building figure 2 with genetic screens data
% SFS and AM 11/18-19/25

%% load datasetVar.mat
load ../'Genetic Screens Dataset'/GeneticScreensDataset.mat
conditionNames = replace(conditionNames,"_", " ");
clear normLFC_noNaN;

% set cutoff Padj
cutoffPadj = 0.25;
%% Removing BHI SMX
% Need to remove BHI condition for SMX so it is balanced since TMP is
% missing condition
tf = contains(conditionNames, 'BHI');
LFC(:,tf) = [];
normLFC(:,tf) = [];
padj(:,tf) = [];
conditionNames(tf) = [];

%% Indexes for SMX and TMP
inxSMX = find(contains(conditionNames, 'SMX'));
inxTMP = find(contains(conditionNames, 'TMP'));

%% Rank LFC
% Calculating rank and percentage LFC
rankLFC = tiedrank(LFC);
pctLFC = round(100 .* rankLFC ./ length(rankLFC));

% Isolate top 5 and bottom 5 percent (or just + padj
logHits = zeros(length(pctLFC), size(pctLFC,2)); % storing logical values
logHitsDep = zeros(length(pctLFC), size(pctLFC,2)); % most sensitive only
for iScreen = 1:size(pctLFC,2) % loop through all screens
    log1 = pctLFC(:,iScreen)<=5 & padj(:,iScreen)<=cutoffPadj; % top 5% most depleted
    log2 = pctLFC(:,iScreen)>=95 & padj(:,iScreen)<=cutoffPadj; % top 5% most enriched
    log = log1 + log2; % combining depleted and enriched
    logHits(:,iScreen) = log; % saving in larger variable
    logHitsDep(:,iScreen) = log1; % saving depleted only
end

% Calculating number of times a gene is a hit
for iGene = 1:length(pctLFC)
    hitsSum(iGene) = sum(logHits(iGene,:));
end

% only want hits that have 2 or more appearances in screens
rankIdx = find(hitsSum>=2);

%% Fig 2B - heat map of correlations between strains and histograms
% Calculate correlation coefficients
[rCond pCond] = corr(pctLFC(rankIdx,:),'rows','pairwise');

% Visualize in heat map
% Set up order by clustering
dCond = pdist(rCond);
treeConds = linkage(dCond,'weighted');
inxCondsCorr = optimalleaforder(treeConds,dCond);

% Plot heat map
fh = figure('color','white');
imagesc(rCond(inxCondsCorr,inxCondsCorr),[-1 1]);
colormap(redbluecmap);
set(gca,'xtick',[1:length(rCond)],'xticklabel',conditionNames(inxCondsCorr));
set(gca,'ytick',[1:length(rCond)],'yticklabel',conditionNames(inxCondsCorr));
axis square; axis tight;
colorbar;
title('correlation between conditions by rank hits');
saveas(fh,'F2 heat map.png');
saveas(fh,'F2 heat map', 'svg');

% Histograms
% Isolate groups to make histograms and remove self correlations
% SMX-SMX
rSMX_SMX = rCond(inxSMX,inxSMX);
rSMX_SMX(rSMX_SMX == 1) = []; % now in one dimension as well

% TMP-TMP
rTMP_TMP = rCond(inxTMP,inxTMP);
rTMP_TMP(rTMP_TMP == 1) = [];  % now in one dimension as well

% SMX-TMP
rSMX_TMP = rCond(inxTMP,inxSMX);

% Plotting histograms
edges = [-.75:.075:.75];
f2 = figure; hold on;
% SMX-SMX histogram
subplot(3,1,1)
histogram(rSMX_SMX,edges, 'Normalization','probability');
xlim([-1 1]);
title('SMX-SMX');
ylim([0 0.35]);
yticks([0:0.05:0.35]);
grid on;

% TMP-TMP histogram
subplot(3,1,2)
histogram(rTMP_TMP,edges, 'Normalization','probability');
xlim([-1 1]);
title('TMP-TMP');
ylim([0 0.35]);
yticks([0:0.05:0.35]);
grid on;

% SMX-TMP histogram
subplot(3,1,3)
histogram(rSMX_TMP,edges, 'Normalization','probability');
xlim([-1 1]);
title('SMX-TMP');
ylim([0 0.35]);
yticks([0:0.05:0.35]);
grid on;

% Save figures
saveas(f2,'F2B hist.png');
saveas(f2,'F2B hist', 'svg');

%% Fig 2C: pie charts figure
% Calculating sums of number of time a k.o. is hit with depleted hits only
c = sum(logHitsDep(rankIdx,:)==1,2); % both SMX and TMP
c(c==0) = []; % removing zeroes to make it truly hits only

% Pie charts all strains
edges = [1 3 5 10 18]-0.1;
h3 = histcounts(c,edges);
fPie2 = figure;
pie(h3)

% Create legend
lgd = legend({'hits in 1-2 screens','hits in 3-4 screens','hits in 5-9 screens', 'hits in 10+ screens'});
colormap gray

% Save figure
saveas(fPie2,'F2C all drugs.png')
saveas(fPie2,'F2C all drugs.svg')

%% Identifying most different genes between SMX and TMP
% Rank sum for comparing SMX and TMP for each gene
p = [];
for i = 1:length(rankIdx) % loop through rank LFC of hits only
    [p(i), h] = ranksum(rankLFC(rankIdx(i), inxSMX), rankLFC(rankIdx(i), inxTMP)); % calculate ranksum p values
end

% False discovery rate calculation
fdr = mafdr(p,'BHFDR',true); % count 243
diffIdx = find(fdr<=0.25); % length 243

% Identify these hits by normLFC instead of rankLFC
rankNormLFC = normLFC(rankIdx,:);
diffNormLFC = rankNormLFC(diffIdx,:);

% Additional cutoff - ensuring different behaviors in each drug
% Calculating median within each index
medsTMP = []; medsSMX = [];
for iGene = 1:length(diffNormLFC)
    medsSMX(iGene) = median(diffNormLFC(iGene,inxSMX));
    medsTMP(iGene) = median(diffNormLFC(iGene,inxTMP));
end

% Checking for 1) gene differs sizably from 0 for one drug and
% 2) other drug has no response or opposite response
count = 1; finIdx = [];
for iGene = 1:length(diffNormLFC)
    medians = [medsSMX(iGene), medsTMP(iGene)];
    % Check that at least one value has norm LFC > 1.5
    if abs(medians(1))> 1.5 || abs(medians(2))>1.5
        % Find less extreme drug response
        [~, absMinIdx] = min(abs(medians)); % identifies SMX or TMP idx as min
        % Check if median is near zero
        if -1 <= medians(absMinIdx) && medians(absMinIdx) <= 1
            finIdx(count) = iGene;
            count = count + 1;
            % Check if median is in opposite direction as max
        elseif medians(1)/medians(2) < 0
            finIdx(count) = iGene;
            count  = count + 1;
        else
        end
    else
    end
end

% Get list of different strains
rankStrains = strains(rankIdx); % first need to get strains by rank idx
diffGenes = rankStrains(diffIdx); % then strains that were significant by rank sum
finDiffGenes = diffGenes(finIdx); % then finally strains that passed last 2 filters

%% Save variables for supp table 1
% Isolating p values to keep
diffFDR = fdr(diffIdx);

save("forTS1", 'diffNormLFC', 'finDiffGenes', 'finIdx', 'inxSMX', 'inxTMP', 'medsSMX', 'medsTMP', 'diffFDR');

%% Fig 2E - Comparing conditions for genes of interest
% Strains of interest
strainsInterest = {'pabA', 'pabB', 'pabC', 'gcvA', 'gcvP', 'gcvT', 'yggS', 'yggR', 'tdk', 'metE', 'pdxA', 'dedD', 'ydgI'};
strainsInterestIdx = [];

% Identify indexes for strains of interest
for iStrain = 1:length(strainsInterest)
    strainsInterestIdx(iStrain) = find(strcmp(strainsInterest{iStrain}, strains));
end

% Plot comparison of SMX and TMP for each gene of interest
% Create new subfolder
if ~exist('F2Eplots/', 'dir')
    mkdir('F2Eplots/');
end

% Plot
for iPlot = 1:length(strainsInterestIdx)
    f2E = figure; hold on;
    SMX = normLFC(strainsInterestIdx(iPlot), inxSMX); % all SMX cond
    TMP = normLFC(strainsInterestIdx(iPlot), inxTMP); % all TMP cond
    tempData = [SMX', TMP'];

    b = boxchart(tempData);

    xlabel('drug');
    xticklabels({'SMX', 'TMP'});
    ylabel('norm LFC');
    title(strainsInterest{iPlot});
    grid on; box on;
    %ylim([0 100]);
    %xlim([0.5 2.5]);
    %xticks([1 2]);

    set(gcf, 'Position', [200 200 450 325]);

    % Save figures
    fileName = fullfile('F2Eplots/', [strainsInterest{iPlot}]);
    saveas(f2E, fileName); % matlab figure
    saveas(f2E, fileName, 'png');
    saveas(f2E, fileName, 'svg');
end
