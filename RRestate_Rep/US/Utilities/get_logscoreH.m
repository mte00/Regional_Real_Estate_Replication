function vecLogScores = get_logscoreH(matPaths, vecData, vecHorizons)
% -------------------------------------------------------------------------
% get_logscore:
%
% Calculates log-scores for a set of predicted values and actual
% observations on a specific variable
% 
% INPUTS:   - matPaths: (# draws)*(# periods) matrix of simulated values for
%               the variable.
%           - vecData: (# periods)*1 vector of observations on the variable
%           - vecHorizons: vector defining the horizons at which the
%               forecasts are evaluated
%
% OUTPUTS:  - vecLogScores: 1*... vector of log-scores. The value
%               is NaN for the densities for which we do not have
%               observations (ie the last ones).
%
% P Alessandri, Jan 2014. Note: this version does not (need to) disciminate
% between "norm" and "cumul"
% -------------------------------------------------------------------------


% Initialise vector of log-scores:
vecLogScores = NaN(1, length(vecHorizons));

% Loop over horizons to get log-score for each predictive density

for i = 1:length(vecHorizons)
    
    ii=vecHorizons(i);
    
    paths = matPaths(:,ii);
    obs   = vecData(ii);
    
    temp=ksdensity(paths,obs);
    if temp==0
        vecLogScores(i)=log(1+temp);
    else
    vecLogScores(i) = log(temp);
    end

end

end