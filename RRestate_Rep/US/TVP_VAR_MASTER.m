% MASTERFILE FOR TVP VAR MODELS TO GET IRFs, FEVDs and network connections
% from REGIONAL MODELS

% MODEL CONTAINS

% REGIONAL RETURNS (MW, NE, S, W) ILLIQ MEASURES (MW, NE, S, W), y, pi, r,
% u ,w, morg_r, hous_st, market_return, mkt_ILLIQ, EPU.

% Sample from 1985M1:2018M12 ALL DATA ARE MONTHLY
addpath('Results')
addpath('Utilities','Data')
clear all; close all; clc;
tic;
% NO NEED TO TRANSFORM DATA ALL IS DONE PRIOR TO ESTIMATING MODELS
for jj=1:4
if jj==1
data=xlsread('US_DATA_VAR_MODELS','R_RtoV_gr','B86:S493');
DFILE='<INSERT PATH HERE>\RRestate_Rep\US\Results\R_RtoV_gr';  
 
elseif jj==2
data=xlsread('US_DATA_VAR_MODELS','N_RtoV_gr','B86:M493');    
DFILE='<INSERT PATH HERE>\RRestate_Rep\US\Results\N_RtoV_gr';

elseif jj==3
data=xlsread('US_DATA_VAR_MODELS','R_IV_gr','B86:S493');  
DFILE='<INSERT PATH HERE>\RRestate_Rep\US\Results\R_IV_gr';  
 
elseif jj==4 % jj=1:4 are models in growth rates
data=xlsread('US_DATA_VAR_MODELS','N_IV_gr','B86:M493');
DFILE='<INSERT PATH HERE>\RRestate_Rep\US\Results\N_IV_gr'; 

end
% ESTIMATE MODEL AND SAVE PARAMETERS AND COVARIANCE MATRICES
Estimate_TVP_VAR_QBLL_singlestep;
%connect; 
jj
clear data X dat2 priorprec0 SI RI a N K L weights w TIC NDC PWC IRF FEV NNN TTT...
    tic ndc pw irfs fevd TT LL NS 
end


