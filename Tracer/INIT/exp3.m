clear all;
modelsetting;
% gives the solution (observations) at the fist step
%CH is the square root o fthe cov matrix, used if generating realizations
%tSols solution
%tObs the observations

% iSol holds the true initial conditions for the true simulation --> true heterogeneity for the FULL domain 

[tSols,tObs] = truemodel(iSol,param);
