function [coeffs,absICarray,fh2, xpoints, ypoints] = calcDoseResponse(drugConc,normAbsorbance,tfPlot,requestedICs, condsUnique, normAbsMG, MGdrugConc, tempFolder)
% function [coeffs,relICarray,absICarray,fh1,fh2, xpoints, ypoints] = calcDoseResponse(drugConc,normAbsorbance,tfPlot,requestedICs)

% Doseresponse function
% Use this function to calculate relative and absolute IC/EC values: 
% Computes the Hill Coefficient and IC50 of a
% dose/response relationship given two vectors describing the doses and
% responses. A semilog graph is plotted illustrating the relationship. The IC50 is also
% labelled. Requires nlinfit from statistics toolbox.
% IC: inhibitory concentrations 
% Absolute:concentration required to bring out the value you want between
% the 0 and maximum response 

% Inputs:
% drugConc: Provide an array of drug concentrations used. Must be horizontal. If you include a
% zero concentration, it will be replaced by four folds below.
% normAbsorbance: Normalized(Blank Subtracted and divided by no drug) Absorbance 600nm values. Must be horizontal.  
% tfPlot: If true a plot will be returned
% requestedICs: Enter a vector of desired IC values (e.g. IC25, IC75 -> [25,75])
% 
% Outputs: 
% coeffs: Coefficients of the fit
% absICarray: Absolute IC Array
% - First row: Desired IC value e.g. IC25->25
% - Second row: IC value concentration
% - Third row: response at IC value 
% fh2: figure handles 
% xpoints, ypoints: values used to plot the sigmoid curve

%% Fitting the sigmoid equation

% Checking if drugConc and norm_abs is vertical and transposing if it is
% not
[r, ~] = size(drugConc);
if r ~=1
    drugConc = drugConc'; 
end 
[r1, ~] = size(normAbsorbance);
if r1 ~=1
    normAbsorbance = normAbsorbance'; 
end 

% Repeat for MG control (or control without EA for EA experiments)
[r, ~] = size(MGdrugConc);
if r ~=1
    MGdrugConc = MGdrugConc'; 
end 
[r1, ~] = size(normAbsMG);
if r1 ~=1
    normAbsMG = normAbsMG'; 
end 

% Hill equation sigmoid (will be used for the dose curve fitting)
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

% Dealing with drug conc 0 (0 conc. will be replaced with very small conc.)
if ~isempty(find(drugConc==0,1))
   drugConcSorted = sort(drugConc,'ascend');
   fc = drugConcSorted(1,3)./drugConcSorted(1,2); % infer the conc jump by 2nd and 3rd conc.
   drugConc(1,drugConc==0) = drugConcSorted(1,2)./(fc^4); % swap with conc. smaller by 4-fold jumps 
end  

% Repeated for control condition (MG or no EA)
% Dealing with drug conc 0 (0 conc. will be replaced with very small conc.)
if ~isempty(find(MGdrugConc==0,1))
   MGdrugConcSorted = sort(MGdrugConc,'ascend');
   fc = MGdrugConcSorted(1,3)./MGdrugConcSorted(1,2); % infer the conc jump by 2nd and 3rd conc.
   MGdrugConc(1,MGdrugConc==0) = MGdrugConcSorted(1,2)./(fc^4); % swap with conc. smaller by 4-fold jumps 
end 

% Calculate some rough guesses for initial parameters
minResponse = min(normAbsorbance);
maxResponse = max(normAbsorbance);
midResponse = mean([minResponse maxResponse]);
minConc = min(drugConc);
maxConc= max(drugConc);

% Repeated for control condition (MG or no EA)
% Calculate some rough guesses for initial parameters
MGminResponse = min(normAbsMG);
MGmaxResponse = max(normAbsMG);
MGmidResponse = mean([MGminResponse MGmaxResponse]);
MGminConc = min(MGdrugConc);
MGmaxConc= max(MGdrugConc);

% Fit the curve and to a sigmoid function 
[coeffs,r,J]=nlinfit(drugConc,normAbsorbance,sigmoid,[minResponse maxResponse midResponse 0]);
IC50=coeffs(3);
hillCoeff=coeffs(4);

% Repeated for control condition (MG or no EA)
[MGcoeffs,MGr,MGJ]=nlinfit(MGdrugConc,normAbsMG,sigmoid,[MGminResponse MGmaxResponse MGmidResponse 0]);
MGIC50=MGcoeffs(3);
MGhillCoeff=MGcoeffs(4);

%% Calculate requested ICs
% Generate sigmoidal fit
NumberXpoints = 5000; % number of points to generate for xpoints 
xpoints=logspace(log10(minConc),log10(maxConc),NumberXpoints); 
ypoints = sigmoid(coeffs,xpoints);

% Repeated for control condition (MG or no EA)
MGxpoints=logspace(log10(MGminConc),log10(MGmaxConc),NumberXpoints); 
MGypoints = sigmoid(MGcoeffs,MGxpoints);

for i=1:length(requestedICs)
        curIC = requestedICs(i); % requested is 50

        % Calculate absolute ICs
        absResponseSpan = max(ypoints)-0; % highest response to drug (most resistant)
        absResponseAtIC =  ((100-curIC)/100)*absResponseSpan+0; % 50% inhibition y axis
        [~,absInx] = min(abs(ypoints-absResponseAtIC)); % index of that 50% inhibition
        absCurDose = round(xpoints(absInx),2); % IC50 calc - x axis value at 50% of y value
        if absInx == NumberXpoints % If parameters not correct for IC50 calc, they = NaN
            absCurDose = nan;
            absResponseAtIC = nan;
        end 
        absICarray(1,i) = curIC; % 50
        absICarray(2,i) = absCurDose; %IC50
        absICarray(3,i) = absResponseAtIC; % y axis

        % Repeated for control condition (MG or no EA)
        absResponseSpanMG = max(MGypoints)-0;
        absResponseAtICMG =  ((100-curIC)/100)*absResponseSpanMG+0;
        [~,absInxMG] = min(abs(MGypoints-absResponseAtICMG));
        absCurDoseMG = round(MGxpoints(absInxMG),2);
        if absInxMG == NumberXpoints
            absCurDoseMG = nan;
            absResponseAtICMG = nan;
        end 
        absICarrayMG(1,i) = curIC;
        absICarrayMG(2,i) = absCurDoseMG;
        absICarrayMG(3,i) = absResponseAtICMG;
end

%% Plot the fitted sigmoid
if(tfPlot) % if instructed to plot
    fh2=figure;
    hold on;
    grid on; 
    % Plot exact points and sigmoid curve
    semilogx(xpoints,ypoints,'Color',[0 0 0],'LineWidth',2); hold on;
    plot(drugConc(1:10), normAbsorbance(1:10), "k.", "LineWidth", 2,"MarkerSize", 40); hold on;
    % Repeated for control condition (MG or no EA)
    semilogx(MGxpoints,MGypoints,'Color',[0 0.5 0.5],'LineWidth',2); hold on;
    plot(MGdrugConc(1:10), normAbsMG(1:10), "b.", "LineWidth", 2, "MarkerSize", 40); hold on;
    yline(1, '--');
    xline(absICarray(2,1));
    xline(absICarrayMG(2,1), 'Color', [0 1 0]);
    xlabel("log(Dose) (uM)")
    ylabel("Normalized Absorbance A600nm")
    ylim([0 1.05]); %ylim([-0.2 1.2])
    xminVar = min(min(drugConc(1:10)), min(MGdrugConc(1:10)));
    xlim([xminVar-(xminVar./2), max([max(xpoints), max(MGxpoints)])]);
    title("Dose Reponse: Absolute IC50");
    % for i=1:length(requestedICs)
    %    curSTR = ['absIC' num2str(absICarray(1,i)) '=' num2str(absICarray(2,i))];
    %     text(absICarray(2,i),absICarray(3,i),[' \leftarrow ' sprintf('%s',curSTR)],'FontSize',12);
    % end
    set(gca, 'xscale', 'log')
    box on;
    saveas(fh2, [tempFolder '/' condsUnique 'wMG.fig']);
    saveas(fh2, [tempFolder '/' condsUnique 'wMG.svg']);
else
    fh2 = 0; 
end 

end
