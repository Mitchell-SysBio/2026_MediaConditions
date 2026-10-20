function [wellCenters, phenotypes, wellNames, h, phenotypesReclass, reClassWell, reClass]=classifyPlatePhenotypes(img, curString, tfRotate,tfFlip)
%function [wellCenters, radius, img2]=classifyPlatePhynotypes(img,strName, tfRotate,tfFlip)
% This function gets an image of a 96-well plate and uses an existing
% neural network output (bestNet) to classify the phenotype of each well

% Setting initial variables
load('bestNet.mat', 'net'); % load the neural network classifier
n_NN_input_layer = 50; % size of input layer for neural network classifier
allPhenotypes = net.Layers(end).Classes; % a list of all phenotypes the classifier knows
myColorMap = cool(length(allPhenotypes)-2);
%myColorMap = [0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1];

myColorMap = [0 0 0; 0 0 0; ... % black (messy, bubbles) cond 1 and 2
    0.9290    0.6940    0.1250; ... % yellow (empty) cond 3
    0.9290    0.6940    0.1250; ... % green (tiny pellet) cond 4
    0.6350    0.0780    0.1840; ... % red (cloudy small) cond 5
    0.6350    0.0780    0.1840; ... % red (cloud white) cond 6
    0.0000    0.4470    0.7410; ... % blue (pellet  small) cond 7
    0.0000    0.4470    0.7410; ... % blue (pellet med) cond 8
    0.6350    0.0780    0.1840; ... % red (cloud large) cond 9
    0.4940    0.1840    0.5560; ... % purple (cloud and pellet large) cond 10
    0.0000    0.4470    0.7410; ... % blue (pellet cloudy) cond 11
    0.0000    0.4470    0.7410;]; % blue (pellet large) cond 12

myAlphaMap = [0.3; 0.3; 0.2; 0.6; 0.2; 0.2; 0.2; 0.2; 0.4; 0.4; 0.4; 0.6];

% Set up 96 well plate dimensions
idRows = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
idColumns = {'1','2','3','4','5','6','7','8','9','10','11','12'};
wellNames = cell(96,1);
wellCenters = nan(96,2);

% Identify whether to rotate and flip the image
if(nargin<4), tfFlip = true; end
if(nargin<3), tfRotate = true; end

% Adjust image (rotate and flip)
if(tfFlip), img = fliplr(img); end
if(tfFlip), img = imrotate(img,90); end

img2 = imadjust(img);
resizeFactor = 1000/min(size(img2)); % adjust the image so the short axis is 1000 pixels
img2 = imresize(img2,resizeFactor);

% Display the image
h = figure('color','w');
movegui(h,"north")
imshow(img2); hold on;
title([curString ': grey = junk; yellow = empty; green = tiny pellet; blue = pellet; purple = cloud/pellet red = clouds']);


% Ask user to select the upper left and bottom right wells
[x, y] = ginput(2); % Get 2 points from the user

% Define the grid dimensions for a 96-well plate (8 rows, 12 columns)
rows = 8;
cols = 12;

% Calculate the distances between wells
rowDist = (y(2) - y(1)) / (rows - 1);
colDist = (x(2) - x(1)) / (cols - 1);

% Draw the grid
radius = round(0.325*mean([rowDist,colDist]));
from_i_to_rc = nan(96,2);
wellNames = cell(96,1);

for row = 0:rows-1
    for col = 0:cols-1
        i = row*12+col+1;
        wellNames{i} = [idRows{row+1} idColumns{col+1}];
        wellCenters(i,1) = x(1) + col * colDist;
        wellCenters(i,2) = y(1) + row * rowDist;
        wellCenters(i,3) = row+1;
        wellCenters(i,4) = col+1;

        from_i_to_rc(i,[1 2]) = [row+1,col+1]; % record conversion from i index to row and col index
        %rect = [wellCenters(i,1)-radius, wellCenters(i,2)-radius, 2*radius, 2*radius];
    end
end

%% Classifing individual wells
% rowSign = ["A", "B", "C", "D", "E","F","G","H"];
% colSign = ["01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"];
l = radius*2;
mask =strel('disk',radius,0);
phenotypes = categorical([]); % place holder for phenotypes

% Resizing images and then classifying
for i = 1:96
    rect = [wellCenters(i,1)-radius, wellCenters(i,2)-radius, l, l];
    curWellImage= imcrop(img2, rect);
    curWellImage(~mask.Neighborhood) = 0;
    img = imresize(curWellImage,[n_NN_input_layer n_NN_input_layer]); % resize image to match proportions of NN classifier
    phenotypes(i) = classify(net, img); % classifying phenotypes using net
end

% Creating image to show user phenotypes for each well
radius = radius*1.25;
for i=1:length(allPhenotypes)
    curPhenotype = allPhenotypes(i);
    tfFitPhenotype = (phenotypes == curPhenotype); % Logical indexing
    % Coloring/marking wells differently based on extent of growth in well
    for j=1:length(tfFitPhenotype)
        if(tfFitPhenotype(j))
            rect = [wellCenters(j,1)-radius, wellCenters(j,2)-radius, 2*radius, 2*radius];
            if(i<=4) % typically considered no growth
                rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1.5,'LineStyle','-');
            elseif(i<=8) % typically considered weak growth
                rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 4,'LineStyle',':');
            else % full growth
                rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 4,'LineStyle','-');
            end

            x = [rect(1), rect(1), rect(1) + rect(3), rect(1) + rect(3)];
            y = [rect(2), rect(2) + rect(4), rect(2) + rect(4), rect(2)];
            fill(x, y,myColorMap(i,:), 'FaceAlpha', myAlphaMap(i));
            text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','k','EdgeColor','k');
        end
    end
end

%% Reclassify based on user input
[reClassWell, reClass] = inputDialog;
phenotypesReclass = phenotypes;
if ~isempty(reClass{1})
    newClass = str2num(reClass{1});
    wellChanged = strsplit(reClassWell{1});
    for reClassed = 1:length(newClass)
        phenotypesReclass(wellNames == string(wellChanged{reClassed})) = allPhenotypes(newClass(reClassed));
    end
    % Changing colormaps to match reclassed inputs
    rectAll = findall(h,'Type', 'Patch');
    set(rectAll, 'Visible', 'off'); % removing patch color
    for i=1:length(allPhenotypes)
        curPhenotype = allPhenotypes(i);
        tfFitPhenotype = (phenotypesReclass == curPhenotype); % Logical indexing
        % Same color labeling method as before, but now updated by user inputs
        for j=1:length(tfFitPhenotype)
            if(tfFitPhenotype(j))
                rect = [wellCenters(j,1)-radius, wellCenters(j,2)-radius, 2*radius, 2*radius];
                if(i<=4) % typically considered no growth
                    rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 1.5,'LineStyle','-');
                elseif(i<=8) % typically considered weak growth
                    rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 4,'LineStyle',':');
                else % full growth
                    rectangle('Position', rect,'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 4,'LineStyle','-');
                end

                x = [rect(1), rect(1), rect(1) + rect(3), rect(1) + rect(3)];
                y = [rect(2), rect(2) + rect(4), rect(2) + rect(4), rect(2)];
                fill(x, y,myColorMap(i,:), 'FaceAlpha', myAlphaMap(i));
                text(x(1),y(1),num2str(i),'BackgroundColor','w','Color','k','EdgeColor','k');
            end
        end
    end
end
%% Saving Quality Assurance Image
saveas(h,['./QAimages/' curString ' - QA.png']);



