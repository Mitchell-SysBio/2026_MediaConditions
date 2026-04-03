function [coeffs,absICarray,fh2, xpoints, ypoints] = calcDoseResponse(drugConc,normAbsorbance,tfPlot,requestedICs, condsUnique, normAbsMG, MGdrugConc, tempFolder)
% function [coeffs,relICarray,absICarray,fh1,fh2, xpoints, ypoints] = calcDoseResponse(drugConc,normAbsorbance,tfPlot,requestedICs)

% doseresponse function
% Use this function to calculate relative and absolute IC/EC values: 
% Computes the Hill Coefficient and IC50 of a
% dose/response relationship given two vectors describing the doses and
% responses. A semilog graph is plotted illustrating the relationship. The IC50 is also
% labelled. Requires nlinfit from statistics toolbox.
% IC: inhibitory concentrations 
% EC: effective concentrations
% Relative: concentration required to bring out the value you
% want between the top and bottom plateaus of the curve 
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
% relICarray: Relative IC Array
% - First row: Desired IC value e.g. IC25->25
% - Second row: IC value concentration
% - Third row: response at IC value 
% absICarray: Absolute IC Array
% - First row: Desired IC value e.g. IC25->25
% - Second row: IC value concentration
% - Third row: response at IC value 
% fh1 and fh2: figure handles 
% xpoints, ypoints: values used to plot the sigmoid curve

%% Fitting the sigmoid equation

% checking if drugConc and norm_abs is vertical and transposing if it is
% not
[r, ~] = size(drugConc);
if r ~=1
    drugConc = drugConc'; 
end 
[r1, ~] = size(normAbsorbance);
if r1 ~=1
    normAbsorbance = normAbsorbance'; 
end 

% NEW repeating for MG ones
[r, ~] = size(MGdrugConc);
if r ~=1
    MGdrugConc = MGdrugConc'; 
end 
[r1, ~] = size(normAbsMG);
if r1 ~=1
    normAbsMG = normAbsMG'; 
end 

% hill equation sigmoid (will be used for the dose curve fitting)
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

% dealing with drug conc 0 (0 conc. will be replaced with very small conc.)
if ~isempty(find(drugConc==0,1))
   drugConcSorted = sort(drugConc,'ascend');
   fc = drugConcSorted(1,3)./drugConcSorted(1,2); % infer the conc jump by 2nd and 3rd conc.
   drugConc(1,drugConc==0) = drugConcSorted(1,2)./(fc^4); % swap with conc. smaller by 4-fold jumps 
end  

% NEW repeated for MG cond
% dealing with drug conc 0 (0 conc. will be replaced with very small conc.)
if ~isempty(find(MGdrugConc==0,1))
   MGdrugConcSorted = sort(MGdrugConc,'ascend');
   fc = MGdrugConcSorted(1,3)./MGdrugConcSorted(1,2); % infer the conc jump by 2nd and 3rd conc.
   MGdrugConc(1,MGdrugConc==0) = MGdrugConcSorted(1,2)./(fc^4); % swap with conc. smaller by 4-fold jumps 
end 

% calculate some rough guesses for initial parameters
minResponse = min(normAbsorbance);
maxResponse = max(normAbsorbance);
midResponse = mean([minResponse maxResponse]);
minConc = min(drugConc);
maxConc= max(drugConc);

% NEW repeated for MG cond
% calculate some rough guesses for initial parameters
MGminResponse = min(normAbsMG);
MGmaxResponse = max(normAbsMG);
MGmidResponse = mean([MGminResponse MGmaxResponse]);
MGminConc = min(MGdrugConc);
MGmaxConc= max(MGdrugConc);

% fit the curve and to a sigmoid function 
[coeffs,r,J]=nlinfit(drugConc,normAbsorbance,sigmoid,[minResponse maxResponse midResponse 0]);
IC50=coeffs(3);
hillCoeff=coeffs(4);

% NEW repeated for MG cond
[MGcoeffs,MGr,MGJ]=nlinfit(MGdrugConc,normAbsMG,sigmoid,[MGminResponse MGmaxResponse MGmidResponse 0]);
MGIC50=MGcoeffs(3);
MGhillCoeff=MGcoeffs(4);

%% calculate requested ICs
NumberXpoints = 5000; % number of points to generate for xpoints 
xpoints=logspace(log10(minConc),log10(maxConc),NumberXpoints); 
ypoints = sigmoid(coeffs,xpoints);

% NEW repeated for MG cond
MGxpoints=logspace(log10(MGminConc),log10(MGmaxConc),NumberXpoints); 
MGypoints = sigmoid(MGcoeffs,MGxpoints);

for i=1:length(requestedICs)
        curIC = requestedICs(i);
        % relResponseSpan = max(ypoints)-min(ypoints);
        % % Calculate relative ICs    
        % relResponseAtIC =  ((100-curIC)/100)*relResponseSpan+min(ypoints); 
        % [~,relInx] = min(abs(ypoints-relResponseAtIC));
        % relCurDose = round(xpoints(relInx),2);
        % if relInx == NumberXpoints
        %     relCurDose = nan;
        %     relResponseAtIC = nan;
        % end 
        % relICarray(1,i) = curIC;
        % relICarray(2,i) = relCurDose;
        % relICarray(3,i) = relResponseAtIC;

        % Calculate absolute ICs
        absResponseSpan = max(ypoints)-0;
        absResponseAtIC =  ((100-curIC)/100)*absResponseSpan+0;
        [~,absInx] = min(abs(ypoints-absResponseAtIC));
        absCurDose = round(xpoints(absInx),2);
        if absInx == NumberXpoints
            absCurDose = nan;
            absResponseAtIC = nan;
        end 
        absICarray(1,i) = curIC;
        absICarray(2,i) = absCurDose;
        absICarray(3,i) = absResponseAtIC;

        % NEW - absolute ICs for MG
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

%% plot the fitted sigmoid

if(tfPlot)
    %anotate the relative IC requested by user
    % fh1=figure('color','white');
    % hold on;
    % grid on; 
    % semilogx(xpoints,ypoints,'Color',[0 0 0],'LineWidth',2);
    % ylim([-0.2 1.2]);
    % plot(drugConc, normAbsorbance, "b*", "LineWidth", 2);
    % xlabel("log(Dose) (uM)");
    % ylabel("Normalized Absorbance A600nm");
    % text(IC50,mean([coeffs(1) coeffs(2)]),[sprintf('rIC_{50}=%0.2f',IC50) ' \rightarrow '],'FontSize',16,'color',[0.5 0.5 0.5],'HorizontalAlignment','right' );
    % title("Dose Reponse: Relative IC50");
    % for i=1:length(requestedICs)
    %     curSTR = ['relIC' num2str(relICarray(1,i)) '=' num2str(relICarray(2,i))];
    %     text(relICarray(2,i),relICarray(3,i),[' \leftarrow ' sprintf('%s',curSTR)],'FontSize',12);
    % end
    %set(gca, 'xscale', 'log')
    %anotate the absolute IC requested by user
    fh2=figure;
    hold on;
    grid on; 
    semilogx(xpoints,ypoints,'Color',[0 0 0],'LineWidth',2); hold on;
    plot(drugConc(1:10), normAbsorbance(1:10), "k.", "LineWidth", 2,"MarkerSize", 40); hold on;
    % MG
    semilogx(MGxpoints,MGypoints,'Color',[0 0.5 0.5],'LineWidth',2); hold on;
    plot(MGdrugConc(1:10), normAbsMG(1:10), "b.", "LineWidth", 2, "MarkerSize", 40); hold on;
    yline(1, '--');
    xline(absICarray(2,1));
    xline(absICarrayMG(2,1), 'Color', [0 1 0]);
    xlabel("log(Dose) (uM)")
    ylabel("Normalized Absorbance A600nm")
    ylim([0 1.05]); %ylim([-0.2 1.2])
    %xminVar = min([min(xpoints(1:10)), min(MGxpoints(1:10))]);
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
