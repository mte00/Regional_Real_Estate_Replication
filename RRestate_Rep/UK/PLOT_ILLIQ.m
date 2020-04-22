addpath('Results')
addpath('Utilities')
addpath('DATA')
addpath('src')
clear all; close all; clc;
yearlab1=(1998:(1/12):2018+(11/12))';
yearlab2=(1985:(1/12):2018+(11/12))';
% LOAD DATA
ukiv=xlsread('UK_DATA_VAR_MODELS','PLOT_MEAS','B2:K253');
usiv=xlsread('UK_DATA_VAR_MODELS','PLOT_MEAS','N2:Q409');

ukiv2=xlsread('UK_DATA_VAR_MODELS','PLOT_MEAS','Y2:AH253');
usiv2=xlsread('UK_DATA_VAR_MODELS','PLOT_MEAS','S2:V409');

temp1=quantile(ukiv,[0.025 0.5 0.975],2);
temp2=quantile(usiv,[0.025 0.5 0.975],2);

temp4=quantile(ukiv2,[0.025 0.5 0.975],2);
temp3=quantile(usiv2,[0.025 0.5 0.975],2);

figure(222)
subplot(2,2,1)
plotx2(yearlab2,[temp3(:,2),temp3(:,1),temp3(:,3)])
axis tight
shadenber()
title('US V^{-1}')

subplot(2,2,2)
plotx2(yearlab1,[temp4(:,2),temp4(:,1),temp4(:,3)])
axis tight
shadenber()
title('UK V^{-1}')

subplot(2,2,3)
plotx3(yearlab2,[temp2(:,2),temp2(:,1),temp2(:,3)])
axis tight
shadenber()
ylabel('%')
xlabel('Time')
title('% deviation from 1-year MA')

subplot(2,2,4)
plotx3(yearlab1,[temp1(:,2),temp1(:,1),temp1(:,3)])
axis tight
shadenber()
ylabel('%')
xlabel('Time')
title('% deviation from 1-year MA')
legend('95% coverage','Median','Location','SouthOutside')
matlab2tikz('ILLIQ.tex')