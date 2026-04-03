function [pwNames, pvalFDR, pwGenes, pwGenesUniv, pwGenesSet]= checkPwEnrichment (geneSet,geneUniverse,pwAnn,pwSizeRange)
% Performing custom GSEA analysis of group of 'hits' compared to defined
% universe without including LFC or counts
% Recieved by MN and adapted by SFS 12-11-25
%% Perform custom Gene Set Enrichment Analysis
% Establish pathwize size range
if (nargin==3) % setting pathway size range if not specified
    pwSizeRange = [3 1000];
end

% Prepare empty variables
pwGenes = {};pwGenesUniv = {};pwGenesSet = {}; % added 02-08-2023
pwSizes = cellfun(@numel,{pwAnn.genes})';
p=[];
k=1;
% Loop through gene ontologies provided in wrapper
for iPw = 1:length(pwAnn)
    if (pwSizes(iPw) > pwSizeRange(1,1) & pwSizes(iPw) < pwSizeRange(1,2)) % if in acceptable pathway size range
        curPw = pwAnn(iPw).genes;
        if (length(intersect(curPw,geneSet))>1) % if ANY overlap between geneSet and genes in current pathway
            if (length(intersect(geneUniverse,geneSet)) < length(geneSet)) % checking if entire gene set is in universe
                disp 'error: gene set not contained in gene universe'
            end
            g1 = geneSet; % group 1 = geneSet
            g2 = intersect(geneUniverse,curPw); % group 2 = genes of current pathway within universe
            N_nono = length(setdiff(geneUniverse,cat(1,g1,g2))); % in gene universe not pathway
            N_yesyes = length(intersect(g1,g2)); % how many genes are in gene set and pathway
            N_yesno = length(setdiff(g1,g2)); % how many genes are in gene set and NOT pathway
            N_noyes = length(setdiff(g2,g1)); % how many genes are in pathway and NOT gene set
            T = [N_nono, N_noyes; N_yesno, N_yesyes];
            [d(k),p(k)] = fishertest(T); % d = decision, p = pvalue
            testedPW(k)=iPw;
            pwGenes{k} = intersect(g1,g2);
            pwGenesUniv{k} = g2;
            pwGenesSet{k} = curPw;
            k=k+1;
        end
    end
end

% Saving outputs
if (length(p)==0)
    pvalFDR = nan;
    pwNames = 'no tests were made';
else
    pvalFDR = mafdr(p,'bhfdr',1)';
    pwNames = {pwAnn(testedPW).names}';
end
end

