%% For compiling max drug concentration from phenotypeImg
% Need to have used "classifyFullPlates..." to classify BMD plates and
% create phenotypeIMG & drugMap variable that contains 
% phenotypeImg{iFile, 1} = name of file 
% phenotypeImg{iFile, 2} = wellCenters= [x, y, row #, column 3]
% phenotypeImg{iFile, 3} = unedited phenotypes
% phenotypeImg{iFile, 4} = edited phenotypes
% phenotypeImg{iFile, 5} = [reClassWell, reClass]; % storing if user needed to edit the phenotype
% phenotypeImg{iFile,6} = edited phenotypes in a 96 well format

% Load phenotypeIMG
load("Phenotyperesults2_12-May-2024.mat")

% Load drug Map
load("drugMap_10-May-2024.mat")
%% Takes edited phenotype and determines maximum drug concentration that has growth
files = phenotypeImg(:,1);
nFiles = size(files,1);

% Setting up variable to contain Max Growth Conc
drugMaxConcV2= zeros(nDrugs, nFiles);
drugCallsV2= zeros(nDrugs, nFiles);

% Looping through all files to calculate Max Growth Conc
for iFile = 1:nFiles
    phenotype96 = phenotypeImg{iFile,6};
    clickMap = cellfun(@(x) str2num(x(1:2)),string(phenotype96)); % making clickmap by making phenotypes into matrix in plate format
    for iDrug = 1:nDrugs % looping through all drugs
        curDrug = drugs(iDrug);
        inxDrug = find(drugMap==curDrug);
        curConc = concMap(inxDrug);
        curConcSort = sort(curConc); % sorting concentration to account for possibility of concentrations not being in order
        curDrugDilutionFactor = curConcSort(2)/curConcSort(1); % calculating the fold conc steps
        curClickValues = clickMap(inxDrug); % order is low to highest concentration of drug
        
        % Making curClickValues into no growth, marginal growth, full growth 
        % Reminder: classes= [01_messy, 02_bubbles, 03_empty, 04_pellet_tiny, 05_cloudy_small, 06_cloud_white,... 
        % 07_pellet_small, 08_pellet_med, 09_cloud_large, 10_cloud_and_pellet_large, 11_pellet_cloudy_large, 12_pellet_large]
        poorGrowthMediaTypes = {'AU', 'RPMI', 'DMEM'}; % for these media, media that have less than phenotype = 8 normally, there is no 'marginal' growth
        imgName = phenotypeImg{iFile,1};
        if contains(imgName, poorGrowthMediaTypes)
            % Poor Growth Media:
            % 1-4 are no growth
            % 5-12 are full growth
            % full growth  = 2; marginal growth = 1; no growth = 0
            curClickValues(curClickValues<=4) = 0; % making anything 4 and less into 0= no growth
            curClickValues(curClickValues==5) = 2;
            curClickValues(curClickValues==6) = 2;
            curClickValues(curClickValues==7) = 2;
            curClickValues(curClickValues==8) = 2;
            curClickValues(curClickValues==9) = 2;
            curClickValues(curClickValues==10) = 2;
            curClickValues(curClickValues==11) = 2;
            curClickValues(curClickValues==12) = 2;
        else
            % Good Media:
            % 1-4 are no growth
            % 5-7 are marginal growth
            % 8-12 are full growth
            % full growth  = 2; marginal growth = 1; no growth = 0
            curClickValues(curClickValues<=4) = 0; % making anything 4 and less into 0= no growth
            curClickValues(curClickValues==5) = 1;
            curClickValues(curClickValues==6) = 1;
            curClickValues(curClickValues==7) = 1;
            curClickValues(curClickValues==8) = 2;
            curClickValues(curClickValues==9) = 2;
            curClickValues(curClickValues==10) = 2;
            curClickValues(curClickValues==11) = 2;
            curClickValues(curClickValues==12) = 2;
        end

        growthInd = find(curClickValues>0); % finding indexes of all that have growth
        

        % Getting maximum growth value by looping through different possibilities
        if ~isempty(growthInd)
            gIndex = 1;
            maxClickValue = curClickValues(growthInd(gIndex));

            % Handing skip well case: called no growth 
            if ~all(curClickValues(growthInd(gIndex):end)) % if there is a no growth well in lower concentrations
                maxClickValue = 0; 
            end

            % Assigning concentration values
            if maxClickValue == 2 % case #1 - full growth was detected
                curMaxConc = curConc(growthInd(gIndex));
                curCall=1;
            elseif maxClickValue == 1 % case #2 = marginal growth was detected in multiple wells
                if length(growthInd(gIndex:end))>1 && (curClickValues(growthInd(gIndex+1)) == 1) % if there are two or more marginal growth, pick the highest concentration with marginal growth
                    margInd = find(curClickValues(growthInd)==1); % finding where marginal growth wells are
                    curMaxConcMarginal = curConc(growthInd(max(margInd))); % getting the highest concentration well with maximum growth to account for trailing marginal growth wells
                    curMaxConcMarginalStepDown = curMaxConcMarginal/curDrugDilutionFactor;
                    curMaxConc = mean([curMaxConcMarginal curMaxConcMarginalStepDown]);
                    curCall=2;
                else % case # 2.5 = only 1 marginal growth well
                    curMaxConcMarginal = curConc(growthInd(gIndex));
                    curMaxConcMarginalStepDown = curMaxConcMarginal/curDrugDilutionFactor;
                    curMaxConc = mean([curMaxConcMarginal curMaxConcMarginalStepDown]);
                    curCall=2.5;
                end
            else % (aka maxClickValue ==0) case #3 = skip well and no other growth
                curMaxConc = NaN;
                curCall=3;
            end

        else % case #4 - no growth reported - will report as next lowest concentration using dilution factor
            curMaxConc = min(curConc)/curDrugDilutionFactor;
            curCall=4;
        end

        % storing concentration with maximum growth
        drugMaxConcV2(iDrug,iFile) = curMaxConc;
        drugCallsV2(iDrug, iFile) = curCall; % call type 
    end
end 

%% Saving variables 
% Making Maximum growth table 
MaxGrowthtableV2 = array2table(drugMaxConcV2);
MaxGrowthtableV2.Properties.VariableNames= strrep(phenotypeImg(:,1),".png","");
MaxGrowthtableV2.Properties.RowNames=string(drugs);

matFN1 = sprintf("CallResultsV2_%s.mat", datetime("today"));
length(drugMaxConcV2)

save(matFN1,"MaxGrowthtableV2", "drugMaxConcV2", "drugCallsV2", "files","drugs")
