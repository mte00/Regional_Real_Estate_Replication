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
data=xlsread('UK_DATA_VAR_MODELS','R_RtoV_gr','B38:AE289');
DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_gr';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_gr_conn'; 
elseif jj==2
data=xlsread('UK_DATA_VAR_MODELS','N_RtoV_gr','B38:M289');    
DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_gr';
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_gr_conn'; 
elseif jj==3
data=xlsread('UK_DATA_VAR_MODELS','R_IV_gr','B38:AE289');  
DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_gr';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_gr_conn'; 
elseif jj==4 % jj=1:4 are models in growth rates
data=xlsread('UK_DATA_VAR_MODELS','N_IV_gr','B38:M289');
DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_gr'; 
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_gr_conn'; 
%elseif jj==5
%data=xlsread('UK_DATA_VAR_MODELS','R_RtoV_lv','B86:S493');
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_lv';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_lv_conn'; 
%elseif jj==6
%data=xlsread('UK_DATA_VAR_MODELS','N_RtoV_lv','B86:M493');
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_lv';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_lv_conn'; 
%elseif jj==7    
%data=xlsread('UK_DATA_VAR_MODELS','R_IV_lv','B86:S493'); 
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_lv'; 
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_lv_conn'; 
%elseif jj==8 % jj=5:8 are models in levels 
%data=xlsread('UK_DATA_VAR_MODELS','N_IV_lv','B86:M493');
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_lv';
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_lv_conn'; 
end
% ESTIMATE MODEL AND SAVE PARAMETERS AND COVARIANCE MATRICES
Estimate_TVP_VAR_QBLL_singlestep;
%connect; 
jj
clear data X dat2 priorprec0 SI RI a N K L weights w TIC NDC PWC IRF FEV NNN TTT...
    tic ndc pw irfs fevd TT LL NS 
end
%%
%clear; clc;

%for jj=1:8
%    tic;
%if jj==1
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_gr';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_gr_conn'; 
%elseif jj==2   
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_gr';
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_gr_conn'; 
%elseif jj==3 
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_gr';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_gr_conn'; 
%elseif jj==4 % jj=1:4 are models in growth rates
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_gr'; 
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_gr_conn'; 
%elseif jj==5
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_lv';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_RtoV_lv_conn'; 
%elseif jj==6
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_lv';  
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_RtoV_lv_conn'; 
%elseif jj==7     
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_lv'; 
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\R_IV_lv_conn'; 
%elseif jj==8 % jj=5:8 are models in levels 
%DFILE='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_lv';
%DFILE2='D:\Dropbox\LRA_Matlab_2020\UK\Results\N_IV_lv_conn'; 
%end
% ESTIMATE MODEL AND SAVE PARAMETERS AND COVARIANCE MATRICES
%load(DFILE);
%connect; 
%toc;
%jj
%end

