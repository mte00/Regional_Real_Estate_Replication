function drule = get_GW_decision(arrDscores, horizons)
% Calculates the "decision rule" suggested by Giacomini-White 2006, p.1558.
%
% Given a vector of cross-model differences in loss functions dL(t) =
% Loss(model 1, t) - Loss(model 2, t) and a set of instruments available in
% t, X(t),  we regress
%
%           dL(t+h) = b'X(t) + eps
% 
% and store the expected discrepancy bhat'X(t). So the rule is "pick mdel 1
% IFF exp discrepancy>0". 
% The instruments are assumed to be X(t) = [1 dL(t)]. We loop over
% horizons, variables and time (t) - the latter because the regression
% must be re-estimated at every t.

% Get dimensions of dL:
[T N H] = size(arrDscores);

drule = NaN(T, N, H);

% Size of smallest (initial) estimation sample:
T0 = 12;

% Loop over horizons, variables, time:
for hh = 1:H
    for ii = 1:N
        for tt = T0+1 : T
            
            h0 = horizons(hh);
            
            % get tt*1 vector of observations available in tt:
            vecDL = arrDscores(1:tt, ii, hh);
            
            % Split into lags and contemporaneous values
            tempX = [ones(tt-h0, 1) vecDL(1:tt-h0)];
            tempY = vecDL(h0+1:end);
            
            % Regress and get expected value E_tt(dL(tt+h0)):
            beta = tempX\tempY;
            tempFitted = beta' * tempX';

            % Store expected value:
            drule(tt, ii, hh) = tempFitted(end);
            
        end     
    end
end


end
