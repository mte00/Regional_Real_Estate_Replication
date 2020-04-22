function pvalues = get_GW_pvalues(arrDscores, type, horizons)
% get_GW_pvalues :
% calculates pvalues for the Giacomini-White (Econometrica, 2006) tests of 
% the null of equal predictive ability. All based on Giacomini's functions 
% (see below).
% 
% INPUTS        - arrDscores: T*N*H array of differences in forecast accuracy
%                   i.e. [L(model 1) - L(model 2)] where L is some loss
%                   function. The dimensions are
%                   (#periods)(#variables)(#horizons)
% 
%                - type: string specifying whether the test is conditional
%                   or unconditional. The conditional test uses as
%                   instrument a lag of the relevant (difference in) scores.
% 
%               - horizons: 1*H vector of horizons 
% 
% OUTPUT        - pvalues: H*N matrix of pvalues.
% 
% NOTES: key ref is G-W, Tests of conditional predictive ability,
% Econometrica (2006). 
% 
% PA, 5 Feb 13.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get dimensions:
[T N H] = size(arrDscores);

pvalues = NaN(H, N);

% Loop over horizons and variables:
for hh = 1:H
    for ii = 1:N
               
        switch type 
            
            case 'uncond'
                lossdiff = squeeze(arrDscores(:, ii, hh));
                instruments = ones(size(lossdiff)) ;
            case 'cond'
                lossdiff = squeeze(arrDscores(2:end, ii, hh));
                instruments = [ones(size(lossdiff)) squeeze(arrDscores(1:end-1, ii, hh))];
        end
                
        [teststat, critval, pval] = CPAtest(lossdiff,instruments,horizons(hh)) ;
        
        pval = round(pval, -3);
        
        pvalues(hh, ii) = pval;
        
    end
end

end
