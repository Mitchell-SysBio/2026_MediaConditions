%% Supplemental Data 4
% SFS 02-14-26
% Categorizing genes and storing in table
%% Load
load forTS4.mat

%% Sort genes
% Table categories defined by most extreme frequency change in a drug
% condition

% Set up empty tables
smxSensTable = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, 'VariableNames', {'Gene Name', 'Median Norm LFC SMX', 'Median Norm LFC TMP', 'FDR adj p val'}); 
smxResTable = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, 'VariableNames', {'Gene Name', 'Median Norm LFC SMX', 'Median Norm LFC TMP', 'FDR adj p val'});
tmpSensTable = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, 'VariableNames', {'Gene Name', 'Median Norm LFC SMX', 'Median Norm LFC TMP', 'FDR adj p val'});
tmpResTable = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, 'VariableNames', {'Gene Name', 'Median Norm LFC SMX', 'Median Norm LFC TMP', 'FDR adj p val'});

% Loop through knockouts and identify what category they belong to
count = 1;
for iGene = finIdx
    medians = [medsSMX(iGene), medsTMP(iGene)];
    % Find most extreme difference
    [~, absMaxIdx] = max(abs(medians));
    absMax = medians(absMaxIdx); % have to get separately to keep sign
    % Sort by index and sign
    if absMaxIdx == 1 % If SMX
        if absMax > 0 % If resistance
            smxResTable = [smxResTable; {finDiffGenes{count}, medians(1), medians(2), diffFDR(iGene)}];
        else
            smxSensTable = [smxSensTable; {finDiffGenes{count}, medians(1), medians(2), diffFDR(iGene)}];
        end
    else % If TMP
        if absMax > 0 % If resistance
            tmpResTable = [tmpResTable; {finDiffGenes{count}, medians(1), medians(2), diffFDR(iGene)}];
        else
            tmpSensTable = [tmpSensTable; {finDiffGenes{count}, medians(1), medians(2), diffFDR(iGene)}];
        end
    end
    count = count + 1;
end

%% Save as excel table with 4 tabs

% Save table to excel file with 4 sheets for 4 categories
filename = 'SuppTables1.xlsx';
sheetNames = {'SMX Res', 'SMX Sens', 'TMP Res', 'TMP Sens'};
writetable(smxResTable, filename, 'Sheet',sheetNames{1});
writetable(smxSensTable, filename, 'Sheet', sheetNames{2});
writetable(tmpResTable, filename, 'Sheet', sheetNames{3});
writetable(tmpSensTable, filename, 'Sheet', sheetNames{4});
