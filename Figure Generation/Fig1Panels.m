%% Figure 1 panel creation
% Note that Figure 1C is separate because it uses a different dataset
% SFS 12-09-25
%% Load datasets
load ../'BMD Dataset'/MGCFig1Variables.mat
load ../'BMD Dataset'/StdPotFig1Variables.mat

%% Figure 1B
f1b = figure; hold on;
for i = 1:length(conditionNames) % looping through each condition
    x = 1:length(drugNames); % all drugs on x axis
    y = MGC_Medians(:,i); % MGCs on y axis

    % Loop through y values to find if there are two medias with the same
    % MGC for the drug to determine if jittering is necessary
    for j = 1:length(y)
        if sum(round(y(j),2)==round(MGC_Medians(j,:),2))>1 && sum(round(y(j),2)==round(MGC_Medians(j,:),2)) <3 % if they're equal to any other media for the drug (>1 because it counts itself)
            x(j) = x(j)+(-0.1 + (rand(1)*0.25)); % adding random space that can be positive or negative
        elseif sum(round(y(j),2)==round(MGC_Medians(j,:),2))>=3 % if there's more than two - so they can be more spaced apart
            x(j) = x(j)+(-0.25 + (rand(1)*0.45));
        else % if not, no jittering
        end

    % Then loop through y values to give them different colored markers
        p = plot(x(j),y(j), 'LineWidth', 1.5, 'MarkerSize', 10, 'Marker', shapeOpt{mediaTypesIdx(i)}); hold on;
        set(p, 'Color', mediaColors(i,:));
        p.LineStyle = "none";
    end
end

% Miscellaneous labels/plot adjustments
xlabel("drugs");
xticks(1:length(drugNames));
xticklabels(drugNames);
ylabel("MGC median");
grid on;
yscale log;
set(f1b, 'Position', [0 0 850 500]);

% Saving figure
savefig(f1b, ['f1b ' char(date) '.fig']);
set(f1b, 'Renderer', 'painters');
saveas(f1b, ['f1b ' char(date) '.svg']);
%exportgraphics(f1b, 'f1b.eps', 'ContentType', 'vector');

% Can make separate figure to generate legend
% Note: ignore the plot itself, just makes it easier to generate legend
% fLegend = figure; hold on;
% x = 1:21;
% y = MGC_Medians(1,:);
% for i = 1:length(y)
%     p2 = plot(x(i), y(i), 'LineWidth', 1.5, 'MarkerSize', 10, 'Marker', shapeOpt{mediaTypesIdx(i)}); hold on;
%         set(p2, 'Color', mediaColors(i,:));
%         p2.LineStyle = "none";
%     legend(conditionNames);
% end
% savefig(fLegend, ['fLegend ' char(date) '.fig']);
% set(fLegend, 'Renderer', 'painters');
% saveas(fLegend, ['fLegend ' char(date) '.svg']);
% exportgraphics(fLegend, 'fLegend.eps', 'ContentType', 'vector');

%% Figure 1D: Clustergram
% Generate clustergram
cgo = clustergram(Z,'ImputeFun', @knnimpute,'Linkage','weighted','RowLabels',drugNames,'ColumnLabels',conditionNames,'DisplayRange',2.5);
cgo.ColumnLabelsRotate = 45; 
cgo.Colormap = redbluecmap;
f1d = figure; plot(cgo);
set(f1d,'position',[295 160 1049 581]);

% NEED TO SAVE CGO AS FIGURE WITH COLOR BAR MANUALLY
exportgraphics(gcf, 'fig1dTry2.eps', 'ContentType', 'vector');
%print(gcf, 'fig1d.svg', '-painters');
% savefig(f1d, ['f1d ' char(date) '.fig']);
% set(f1d, 'Renderer', 'painters');
% saveas(f1d, ['f1d ' char(date) '.svg']);
% exportgraphics(f1d, 'f1d.eps', 'ContentType', 'vector');

%% Figures 1E and 1F

% Calculate correlations between media and between drugs
[rDrugs pDrugs] = corr(Z','rows','pairwise');
[rMedia pMedia] = corr(Z,'rows','pairwise');

% Set up order of heat maps based on clustering
dDrugs = pdist(rDrugs);
dMedia = pdist(rMedia);
treeDrugs = linkage(dDrugs,'weighted');
treeMedia = linkage(dMedia,'weighted');
inxDrugs = optimalleaforder(treeDrugs,dDrugs);
inxMedia = optimalleaforder(treeMedia,dMedia);

% Produce heat map 1: antibiotics
f1e = figure('color','white');
imagesc(rDrugs(inxDrugs,inxDrugs),[-1 1]); 
colormap(redbluecmap);
set(gca,'xtick',[1:length(rDrugs)],'xticklabel',drugNames(inxDrugs));
set(gca,'ytick',[1:length(rDrugs)],'yticklabel',drugNames(inxDrugs));
axis square; axis tight;
colorbar;

% Saving figure
savefig(f1e, ['f1e ' char(date) '.fig']);
set(f1e, 'Renderer', 'painters');
saveas(f1e, ['f1e ' char(date) '.svg']);
exportgraphics(f1e, 'f1e.eps', 'ContentType', 'vector');

% Produce heat map 2: media
f1f = figure('color','white');
imagesc(rMedia(inxMedia,inxMedia),[-1 1]); 
colormap(redbluecmap);
set(gca,'xtick',[1:length(rMedia)],'xticklabel',conditionNames(inxMedia));
set(gca,'ytick',[1:length(rMedia)],'yticklabel',conditionNames(inxMedia));
axis square; axis tight;
colorbar;

% Saving figure
%set(gcf,'position',[295 160 1049 581]);
savefig(f1f, ['f1f ' char(date) '.fig']);
set(f1f, 'Renderer', 'painters');
saveas(f1f, ['f1f ' char(date) '.svg']);
exportgraphics(f1f, 'f1f.eps', 'ContentType', 'vector');