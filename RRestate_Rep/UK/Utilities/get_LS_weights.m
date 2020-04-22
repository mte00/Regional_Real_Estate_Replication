function matWeights = get_LS_weights(matData1, matData2, strWeightType)
% get_LS_weights
%
% Calculates weights to be used to re-weight a (set of) predictive log
% scores. The weight for time t is only a function of the observation Y(t)
% and the sample on which the model has been estimated [Y(t0), ...,
% Y(t-1)], so the same weights are used for all models.
% 
% Reference: Amisano and Giacomini, Comparing density forecasts via
% weighted likelihood ratio tests, JBES 2007 (AG).
%
% INPUTS:   - matData1: T1*N matrix of observations on N variables,
%               to be used to estimate their uncoditional distribution as
%               of t=T1
%
%           - matData2: H*N matrix of observations on N variables, to be
%               weighted. 
%
%           - strWeightType: type of weights to be calculated, i.e. region 
%               of the density on which the evaluation should focus (see AG
%               for details). Acceptable values: 'centre', 'tails',
%               'ltail', 'rtail'.
%
% OUTPUTS:  - matWeights: H*N matrix of weights for the log-scores attached
%               to the observations matData2
%
% This version: P Alessandri, Jan 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions:
[T N] = size(matData1);
[H N2] = size(matData2);

% if N1 neq N2, error

% Output 
matWeights = NaN(N, H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get standardised data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matMean = repmat(mean(matData1,1), [H 1]);
matStdv = repmat(std(matData1, 0, 1), [H 1]);

matDataStar = ( matData2 - matMean ) ./ matStdv;
% This is Yst in Amisano-Giacomini equation 1, p.179.
% Note that we standardise matData2 with moments calculated using matData1.
% Data1 is supposed to contain lags of the Data2.


%% Calculate weights using standard normal distribution:
switch strWeightType
    
    case 'centre'
        matWeights = pdf('Normal', matDataStar, 0, 1);
        
    case 'tails'
        matWeights = 1 - pdf('Normal', matDataStar, 0, 1) / pdf('Normal', 0, 0, 1) ;
    
    case 'ltail'
        matWeights = 1 - cdf('Normal', matDataStar, 0, 1);
        
    case 'rtail'
        matWeights = cdf('Normal', matDataStar, 0, 1);
        
    otherwise
        error('Unknown weighting function: check inputs')
end


end
