% DENSITY FORECASTING MASTER FILE USING 4-VARIABLE VAR OF MONTHLY US
% HOUSING RETURNS AND MEASURES OF ILLIQUIDITY AS PROPOSED IN EFZ (202X)

clear all; close all; clc; warning off;

strFolder = '<INSERT PATH HERE>\RRestate_Rep\UK';

addpath('Results')
addpath('Data')
addpath('Utilities')
addpath('Forecasts')

% REGIONAL
% LOAD DATA
RETS=xlsread('UK_DATA_VAR_MODELS','R_RtoV_gr','B38:AE289'); 
RETS=[RETS(:,1:10), RETS(:,21:end)];
IVOL=xlsread('UK_DATA_VAR_MODELS','R_IV_gr','B38:AE289'); 
RTOV=xlsread('UK_DATA_VAR_MODELS','R_RtoV_gr','B38:AE289');  
nsim=5000; shrinkage=0.05; L=2;
T_thres=60; % corresponds to rolling window of 48 month horizon
nfore=12;
%
% FIRST DO TVP VAR WITH NO ILLIQ
data=RETS;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'Regional',filesep,'NO_LIQ',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH RtoV
data=RTOV;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'Regional',filesep,'RtoV',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH iV
data=IVOL;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'Regional',filesep,'invV',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);


 Evaluate_BVAR_forecasts;
 
  clear scoreNOLIQ scoreRTOV scoreIVOL pitsNOLIQ pitsRTOV pitsIVOL...
   wlsnoliq_L wlsnoliq_LR wlsrtov_L wlsrtov_LR wlsinvv_L wlsinvv_LR RETS RTOV IVOL

%
%
% NATIONAL
% LOAD DATA
RETS=xlsread('UK_DATA_VAR_MODELS','N_RtoV_gr','B38:M289'); 
RETS=[RETS(:,1), RETS(:,3:end)];
IVOL=xlsread('UK_DATA_VAR_MODELS','N_IV_gr','B38:M289'); 
RTOV=xlsread('UK_DATA_VAR_MODELS','N_RtoV_gr','B38:M289');
nsim=5000; shrinkage=0.05; L=2;
T_thres=60; % corresponds to 60 month initial horizon
nfore=12;
%
% FIRST DO TVP VAR WITH NO ILLIQ
data=RETS;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'National',filesep,'NO_LIQ',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH RtoV
data=RTOV;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'National',filesep,'RtoV',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH iV
data=IVOL;
fdirname=strcat(strFolder,filesep,'BVAR',filesep,'National',filesep,'invV',filesep); %directory to save forecast densities
out=BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

 Evaluate_BVAR_forecasts_2;