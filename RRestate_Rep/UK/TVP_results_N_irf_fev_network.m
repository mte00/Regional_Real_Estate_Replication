clear; clc;
addpath('src')
addpath('Results')
addpath('Utilities')

HO=21; % horizon we want to look at for FEVDS
Nregion=1; % 1 regions for national 
% first do RtoV regional
yearlab=(1998.00:(1/12):2018.75)';

load N_RtoV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);


FE1=zeros(TTT,3);
FEE1=zeros(TTT,12);
parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);

   ir1(kk,:,:)=squeeze(temp1);

   
   fe1=zeros(1,Nregion,3);
   fee=zeros(1,12);
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,1,:,:);    
   fe1(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:12
      temp5=FEV{kk,1}(jj,1,2,:); 
      fee(1,jj)=temp5;
   end
   
   
   FE1(kk,:,:)=fe1;
   FEE1(kk,:)=fee;
end

% IRFs of return for region j wrt to illiq shock region j

D1=120:138; 

IR1_1=squeeze(mean(ir1(D1,:,:),1));



figure(1)
subplot(2,1,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
title('Average Response of UK Return over Great Recession')
ylabel('%')
xlabel('Horizon')
subplot(2,1,2)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
matlab2tikz('UK_IRF_RTOV_N_20.tex')

figure(2)
area(yearlab, 100*FE1(:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('UK')
matlab2tikz('UK_FEVD_ILLIQ_RTOV_N.tex')
%%
% Now for network stuff

% total connectedness
figure(11)
plot(yearlab,TIC(:,2),'b-','LineWidth',1.1)
hold on,
plot(yearlab,TIC(:,1),'b--')
plot(yearlab,TIC(:,3),'b--')
shadenber()
axis tight
xlabel('Time')
matlab2tikz('UK_TIC_RtoV_N.tex')

% net directional regions returns
figure(12)
subplot(2,1,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-3 15])
shadenber()
title('Return')
subplot(2,1,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-3 15])
shadenber()
title('Illiquidity')
matlab2tikz('UK_NDC_RTV_N.tex')





%%

% REGIONAL IVOL NOW
load N_IV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);


FE1=zeros(TTT,3);
FEE1=zeros(TTT,12);
parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);

   ir1(kk,:,:)=squeeze(temp1);

   
   fe1=zeros(1,Nregion,3);
   fee=zeros(1,12);
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,1,:,:);    
   fe1(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:12
      temp5=FEV{kk,1}(jj,1,2,:); 
      fee(1,jj)=temp5;
   end
   
   
   FE1(kk,:,:)=fe1;
   FEE1(kk,:)=fee;
end

% IRFs of return for region j wrt to illiq shock region j

D1=120:138; 

IR1_1=squeeze(mean(ir1(D1,:,:),1));



figure(3)
subplot(2,1,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
title('Average Response of UK Return over Great Recession')
ylabel('%')
xlabel('Horizon')
subplot(2,1,2)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
matlab2tikz('UK_IRF_IV_N_20.tex')

figure(4)
area(yearlab, 100*FE1(:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('UK')
matlab2tikz('UK_FEVD_ILLIQ_IV_N.tex')

%%
% Now for network stuff

% total connectedness
figure(21)
plot(yearlab,TIC(:,2),'b-','LineWidth',1.1)
hold on,
plot(yearlab,TIC(:,1),'b--')
plot(yearlab,TIC(:,3),'b--')
shadenber()
axis tight
xlabel('Time')
matlab2tikz('UK_TIC_IV_N.tex')

% net directional regions returns
figure(22)
subplot(2,1,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-3 15])
shadenber()
title('Return')
subplot(2,1,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-3 15])
shadenber()
title('Illiquidity')
matlab2tikz('UK_NDC_IV_N.tex')
