% This wrapper load images of 96-well plates and uses
% classifyPlatePhynotypes to classify phenotype based on neural network 

% load plate map metadata and identify plate images for analysis
platemapFileName = 'drugConcMetadata.xlsx';
drugMap = categorical(readcell(platemapFileName,'Range','B12:M19'));
concMap = readmatrix(platemapFileName,'Range','B2:M9');

dataFolder = './BMDimages/';
dataFormat = '.png';
files = dir([dataFolder '*' dataFormat]); 
nFiles = size(files,1);

% Prepare place holder for drug sensitivity data
drugs = unique(drugMap); nDrugs = numel(drugs);

%% Classifying wells
% classes= [01_messy, 02_bubbles, 03_empty, 04_pellet_tiny, 05_cloudy_small, 06_cloud_white, 07_pellet_small, 08_pellet_med, 09_cloud_large, 10_cloud_and_pellet_large, 11_pellet_cloudy_large, 12_pellet_large]

% Establish variable for storing the phenotypes
if exist('phenotypeImg') == 0
    intialize phenotypeImg 
   phenotypeImg = cell(nFiles,3); % to store the phenotype 
    phenotypeImg{nFiles,5} =  [];
    spacer = 0; 
    startNum = 1;
else
     spacer = 0; % use if not starting from file = 1 
     startNum = spacer+1;
end 

% Loop through plate images and classify using neural network output
% bestNet
for iFile = 1:(nFiles+spacer)
    curFile = [dataFolder files(iFile-spacer).name];
    tokens = split(curFile,'/');
    curString = tokens{3};
    img = imread(curFile);
    % classify using neural network 
    [wellCenters, phenotypes, wellNames, h, phenotypeReclass, reClassWell, reClass]=classifyPlatePhenotypes(img, curString);
    close(h)
    % store everything in phenotypeImg variable
    phenotypeImg{iFile, 1} = curString; 
    phenotypeImg{iFile, 2} = wellCenters;  
    phenotypeImg{iFile, 3} = phenotypes; % unedited phenotypes
    phenotypeImg{iFile, 4} = phenotypeReclass; % edited phenotypes
    phenotypeImg{iFile, 5} = [reClassWell, reClass]; % storing if user needed to edit the phenotype
    phenotype96 = reshape(phenotypeReclass, 12,8)';
    phenotypeImg{iFile, 6} = phenotype96; % storing phenotypes in 96well format
end 

%% Saving variables 
matFN = sprintf("Phenotyperesults2_%s.mat", datetime("today"));
save(matFN, "phenotypeImg")

% saving drug metadata
matFN2 = sprintf("drugMap_%s.mat", datetime("today"));
save(matFN2, "drugMap", "concMap", "drugs", "nDrugs");
