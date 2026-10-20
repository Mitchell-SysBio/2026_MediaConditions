%% Growth measurement from non-shaking culture
% Initial analysis from experiment to evaluate endpoint growth measurement
% by absorbance of non-shaking culture

%% Bring in data
background = table2array(readtable('background_062024.xlsx')); % blank for each media
endPointGrowth = table2array(readtable('endpointMeasure_062024.xlsx')); % 24 hr time point measurements

media_order = readcell("metadata_062024.xlsx"); % order of media

%% Average and subtract background
% Take average
endPointGrowth_av = mean(endPointGrowth, 'omitnan');

endPt_norm = endPointGrowth_av - background; % subtract blank for each media

% Saving for comparison
save('StatGrowth.mat', 'endPt_norm', 'media_order');