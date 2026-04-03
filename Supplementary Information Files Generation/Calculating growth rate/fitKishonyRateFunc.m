function [sol,f,growthRate, doublingTime, yield, lagTime5perc] = fitKishonyRateFunc(x,y)
%% Fit for bacterial growth curves
% Based on the fitting function from Hegreness et al. PNAS 2008
% "We estimated the growth rate, r [as well as the three parameters for the baseline, a, yield, a + b, and lag time, ln(c)/r]
% for each population on each day by the best fit of the function OD(t) = a + b/[1 + cexp(?r?t)]

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0,0],...
    'Upper',[0.001,inf,inf,inf],...
    'Startpoint',[0 max(y) 1 1]);

f = fittype('a+b/(1+c*exp(-r*x))','options',s); % set fit equation

sol = fit(x',y',f,'Robust','on'); % calculate fit

% Isolate variables from fit function
% a=baseline; a+b=yield; lag time= ln(c)/r; gr = r;
growthRate = sol.r;
yield = sol.a+sol.b;
lagTime5perc= log(19/sol.c)/-sol.r; % originally by calculated as time to 50% of yield (lag time= ln(c)/r)
% we calculate time to 5% of yield (ln(19/c)/-r)
doublingTime = log(2)/sol.r;