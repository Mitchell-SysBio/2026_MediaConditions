%% Fig 3 - Comparing results from GSEA using REVIGO
% SFS 12-19-25
% Takes in results from GSEA analysis to be simplified by Revigo website
% Then reorganizies Revigo output to be used for figures

% Note - NEED TO MANUALLY USE INPUTS FOR REVIGO before completing second
% half of script

%% Load GSEA results
fileName = 'gseaResults.mat';
load(fileName)
annotationID = 1; % go BP

% Assign separate variables for chosen annotation
pwNames= gseaOutput.PwName{annotationID};
genesPresent = gseaOutput.genesPresent{annotationID};
genesUniverse = gseaOutput.genesUniverse{annotationID};
padj = gseaOutput.pvalFDR{annotationID};

%% Preparing for REVIGO
revigoInput = {};
% Looping through all enriched/depleted pathways
for iPW = 1:length(pwNames)
    curEnrichedTemp = strsplit(pwNames{iPW});
    revigoInput{iPW,1} = curEnrichedTemp{1}; % Add pathway identifier
    revigoInput{iPW,2} = padj(iPW); % Add p adj value
end
% Write input to excel document that can be copy and pasted into REVIGO
writecell(revigoInput,'forRevigo.xlsx');

%% Post REVIGO - only run AFTER running revigo manually
% First need to manually compose a variable with revigo results!
% Load manually reordered by value results for chosen conditions
resP01 = readtable('Revigo_PPoint01.xlsx'); % PUT YOUR VARIABLE NAME HERE WITH REVIGO OUTPUT

% Find revigo results in pw names
for iPW = 1:height(resP01) % for all pathways from Revigo
    temp = find(contains(pwNames, resP01.TermID(iPW)));
    if isempty(temp)
        pPWIdx(iPW) = NaN;
    else
        pPWIdx(iPW) = temp; % find index of each Revigo pathway in pwNames variable
    end
end

% Remove results with large universe sizes (>50) due to low specificity
% Isolate universe sizes (revigo results only)
for iPW = 1:length(pPWIdx) % For all revigo pathways found in pwNames
    if ~isnan(pPWIdx(iPW)) % if not empty -> count size of universe
        pUnivCount(iPW) = length(genesUniverse{pPWIdx(iPW)});
    else % if empty -> NaN
        pUnivCount(iPW) = NaN;
    end
end

% Find universe sizes greater than 50
pRemoveIdx = find(pUnivCount>50);

% Remove revigo results with greater than 50 universe size
resP01(pRemoveIdx,:) = []; pUnivCount(pRemoveIdx) = [];

% isolate top 10 revigo hits with under 50 universe size
resP01_10 = resP01(1:10,:);

% Reorder data by max gene universe size
[~, pOrder] = sort(pUnivCount(1:10), 'descend');
pUnivCount_Sorted = pUnivCount(pOrder);
resP01_10_Sorted = resP01_10(pOrder,:);

%% Save variables
save('pwFinVars.mat', 'resP01_10_Sorted', 'pUnivCount_Sorted');
