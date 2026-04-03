%% Pathway enrichment analysis on manually curated datasets
% modified 12-11-25 by SFS, received from MN
%% Load gene variables
geneSet = table2cell(readtable('pCorr.txt',ReadVariableNames=false));
geneUniverse = table2cell(readtable('Universe.txt', ReadVariableNames=false));

% Load pw annottions from Creating Pathway Annotations folder
load 'pwAnnVar.mat';

%% Define cutoffs
pvalFDR_cutoff_strain = 0.1; % cutoff to define a significant enriched strains (from previous analysis)
pvalFDR_cutoff_pw = 0.25; % cutoff to define significantly enriched pathways

%% Other definitions
pwSizeRange = [2 1500];
ontologyArray = {goBpAnn;goCcAnn;goMfAnn;keggAnn;cogAnn;ecocycAnn};
ontologyDesc = {'goBpAnn'; 'goCcAnn'; 'goMfAnn'; 'keggAnn';'cogAnn';'ecocycAnn'};

%% Using PwEnrichment function
for j = 1:length(ontologyArray) % looping through different ontologies
    if (j==6 | j==4)
        pwSizeRange = [2 50]; % change range of hits depending on database characteristics
    end
    [PwName, pvalFDR, PwGenes, PwGenesUniv, PwGeneSet] = checkPwEnrichment(geneSet,geneUniverse,ontologyArray{j},pwSizeRange);
    inx = pvalFDR<pvalFDR_cutoff_pw; % checking that variables match expected pvalFDR significance
    PwName = PwName(inx); pvalFDR = pvalFDR(inx); % overwrite variables keeping only significant ones
    PwGenes = PwGenes(inx); PwGenesUniv = PwGenesUniv(inx); PwGeneSet = PwGeneSet(inx); % continuing with additional variables
    [~,inx] = sort(pvalFDR,'ascend'); % sorting by lowest padj value
    % Assign temporary variables into permanent variable
    gseaOutput.pvalFDR{j} = pvalFDR(inx);
    gseaOutput.PwName{j} = PwName(inx);
    gseaOutput.OntologyName(j) = ontologyDesc(j);
    gseaOutput.genesPresent{j} = PwGenes(inx)';
    gseaOutput.genesUniverse{j} = PwGenesUniv(inx)';
    gseaOutput.genesSet{j} = PwGeneSet(inx)';
end

%% Save gseaOutput
save('gseaResults.mat', 'gseaOutput');
