function [teststat,critval,pval]=CPAtest(lossdiff,instruments,tau)

% This function performs the asymptotic Conditional Predictive Ability Test
% INPUTS: lossdiff, Tx1 vector of differences in losses over the out of sample period for the two models under consideration
%         instruments, a Txk matrix corresponding to the test function ht
%         tau, the forecast horizon
%
% OUTPUTS: teststat, the test statistic of the conditional predictive ability test
%          critval, the critical value of the test for a 5% level (the test is a chi square)
%          pval, the p-value of the test
%
% Raffaella Giacomini, 2003. 
% Source: http://www.homepages.ucl.ac.uk/~uctprgi/ 


T = size(lossdiff,1);
% create the regressor matrix given by lossdiff*ht', where ht is the matrix of instruments
reg = -999*ones(size(instruments));
for jj = 1:size(instruments,2)
   reg(:,jj) = instruments(:,jj).*lossdiff;
end

if tau == 1
   % calculate the test stat as nR^2 from the regression of one on lossdiff*ht
   res.beta = reg\ones(T,1);   
   err = ones(T,1)-reg*res.beta;
   r2 = 1-mean(err.^2);
   teststat = T*r2;
   q = size(reg,2);
   critval = chi2inv(.95,q);
   pval = 1 - cdf('chi2',abs(teststat),q);
else
   zbar = mean(reg)';
   nlags = tau-1;
   omega = NeweyWest_matrix(reg,nlags);
   teststat = T*zbar'*inv(omega)*zbar;
   q = size(reg,2);
   critval = chi2inv(.95,q);
   pval = 1 - cdf('chi2',abs(teststat),q);
end
   
   