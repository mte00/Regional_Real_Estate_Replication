% DENSITY FORECASTING MASTER FILE USING 4-VARIABLE VAR OF MONTHLY US
% HOUSING RETURNS AND MEASURES OF ILLIQUIDITY AS PROPOSED IN EFZ (202X)

clear all; close all; clc; warning off;

strFolder = '<INSERT PATH HERE>\RRestate_Rep\UK';
%strFolder = '/Users/mte00/Dropbox/RRestate_Rep/UK';
% JAN 1998--DEC2018
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
T_thres=60; % corresponds to 60 month initial horizon
nfore=12;
%
% FIRST DO TVP VAR WITH NO ILLIQ
data=RETS;
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'Regional',filesep,'NO_LIQ',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH RtoV
data=RTOV;
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'Regional',filesep,'RtoV',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH iV
data=IVOL;
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'Regional',filesep,'invV',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);
%
 Evaluate_forecasts;
 clear scoreNOLIQ scoreRTOV scoreIVOL pitsNOLIQ pitsRTOV pitsIVOL...
    wlsnoliq_L wlsnoliq_LR wlsrtov_L wlsrtov_LR wlsinvv_L wlsinvv_LR RETS RTOV IVOL
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
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'National',filesep,'NO_LIQ',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH RtoV
data=RTOV;
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'National',filesep,'RtoV',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);

% NEXT DO TVP VAR WITH iV
data=IVOL;
fdirname=strcat(strFolder,filesep,'Forecasts',filesep,'National',filesep,'invV',filesep); %directory to save forecast densities
out=TVPVAR_QBLL_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname);
%
 Evaluate_forecasts_2;
 
 
%
%%
addpath('src');

%%
Nregion=10;
hhp=4;
yearlab=(1998.00:(1/12):2018+(11/12))';

yearlab=yearlab(T0+L:end);
yearlab=yearlab(1:T2,:);

label={'EE';'EM';'NE';'NW';'LO';'SE';'SW';'WA';'WM';'YO'};

for vv=1:Nregion
       % Get marginal likelihoods: ML(t) = SUM_t(LS(t), Geweke-Amisano equ. (8)
    ml1 = []; ml2 = []; ml3 = [];
    for tt = 1:T2
        ml1 = [ml1 ; sum(NO_LIQ.ls(1:tt,vv,hhp)) ];
        ml2 = [ml2 ; sum(inV_LIQ.ls(1:tt,vv,hhp)) ];
        ml3 = [ml3 ; sum(RtV_LIQ.ls(1:tt,vv,hhp)) ];   
    end
    % Get log Bayes factors (Geweke-Amisano eq. 9):
    bf21 = ml2 - ml1;
    bf31 = ml3 - ml1;   
    figure(332)
subplot(5, 2, vv)
    hold on 
    plot(yearlab, bf31, 'b', 'LineWidth', 2)
    plot(yearlab, bf21, 'r--', 'LineWidth', 2)
    plot(yearlab,zeros(1,length(yearlab)),'k-','LineWidth',1.2)
    axis tight
    shadenber()
    title(label(vv));
    if vv==Nregion
    matlab2tikz('BF_UK.tex')
    end
    hold off
end

vv=1;
       % Get marginal likelihoods: ML(t) = SUM_t(LS(t), Geweke-Amisano equ. (8)
    ml1 = []; ml2 = []; ml3 = [];
    for tt = 1:T2
        ml1 = [ml1 ; sum(NO_LIQ_N.ls(1:tt,vv,hhp)) ];
        ml2 = [ml2 ; sum(inV_LIQ_N.ls(1:tt,vv,hhp)) ];
        ml3 = [ml3 ; sum(RtV_LIQ_N.ls(1:tt,vv,hhp)) ];   
    end
    % Get log Bayes factors (Geweke-Amisano eq. 9):
    bf21 = ml2 - ml1;
    bf31 = ml3 - ml1;   
    figure(333)
    hold on 
    plot(yearlab, bf31, 'b', 'LineWidth', 2)
    plot(yearlab, bf21, 'r--', 'LineWidth', 2)
    plot(yearlab,zeros(1,length(yearlab)),'k--','LineWidth',1.2)
    shadenber()
    title('UK')
    legend('RtoV vs No Liq','V^{-1} vs No Liq','Location','SouthOutside')
    hold off
    matlab2tikz('BF_UK_N.tex')
    
    
    
    for vv=1:Nregion
       % Get marginal likelihoods: ML(t) = SUM_t(LS(t), Geweke-Amisano equ. (8)
    ml1 = []; ml2 = []; ml3 = [];
    for tt = 1:T2
        ml1 = [ml1 ; sum(NO_LIQ.wlsL(1:tt,vv,hhp)) ];
        ml2 = [ml2 ; sum(inV_LIQ.wlsL(1:tt,vv,hhp)) ];
        ml3 = [ml3 ; sum(RtV_LIQ.wlsL(1:tt,vv,hhp)) ];   
    end
    % Get log Bayes factors (Geweke-Amisano eq. 9):
    bf21 = ml2 - ml1;
    bf31 = ml3 - ml1;   
    figure(442)
subplot(5, 2, vv)
    hold on 
    plot(yearlab, bf31, 'b', 'LineWidth', 2)
    plot(yearlab, bf21, 'r--', 'LineWidth', 2)
    plot(yearlab,zeros(1,length(yearlab)),'k-','LineWidth',1.2)
    axis tight
    shadenber()
    title(label(vv));
    if vv==Nregion
    matlab2tikz('BF_UK_wls.tex')
    end
    hold off
end

vv=1;
       % Get marginal likelihoods: ML(t) = SUM_t(LS(t), Geweke-Amisano equ. (8)
    ml1 = []; ml2 = []; ml3 = [];
    for tt = 1:T2
        ml1 = [ml1 ; sum(NO_LIQ_N.wlsL(1:tt,vv,hhp)) ];
        ml2 = [ml2 ; sum(inV_LIQ_N.wlsL(1:tt,vv,hhp)) ];
        ml3 = [ml3 ; sum(RtV_LIQ_N.wlsL(1:tt,vv,hhp)) ];   
    end
    % Get log Bayes factors (Geweke-Amisano eq. 9):
    bf21 = ml2 - ml1;
    bf31 = ml3 - ml1;   
    figure(443)
    hold on 
    plot(yearlab, bf31, 'b', 'LineWidth', 2)
    plot(yearlab, bf21, 'r--', 'LineWidth', 2)
    plot(yearlab,zeros(1,length(yearlab)),'k--','LineWidth',1.2)
    shadenber()
    title('UK')
    legend('RtoV vs No Liq','V^{-1} vs No Liq','Location','SouthOutside')
    hold off
    matlab2tikz('BF_UK_N_wls.tex')



