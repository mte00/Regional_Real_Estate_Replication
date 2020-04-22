clear; clc;
addpath('src')
addpath('Results')
addpath('Utilities')
% ORDER OF REGIONS:
% EE EM NE NW LO SE SW WA WM YO

HO=21; % horizon we want to look at for FEVDS
Nregion=10; % 10 regions for UK in our study
% first do RtoV regional
yearlab=(1998.00:(1/12):2018.75)';

load R_RtoV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);
ir2=ir1; ir3=ir1; ir4=ir1; ir5=ir1; ir6=ir1; ir7=ir1; ir8=ir1; ir9=ir1; ir10=ir1;

FE1=zeros(TTT,Nregion,3);
FE2=FE1; FE3=FE1; FE4=FE1; FE5=FE1; FE6=FE1; FE7=FE1; FE8=FE1; FE9=FE1; FE10=FE1;

parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);
   temp2=IRF{kk,1}(Nregion+2,2,:,:);
   temp3=IRF{kk,1}(Nregion+3,3,:,:);
   temp4=IRF{kk,1}(Nregion+4,4,:,:);
   temp5=IRF{kk,1}(Nregion+5,5,:,:);
   temp6=IRF{kk,1}(Nregion+6,6,:,:);
   temp7=IRF{kk,1}(Nregion+7,7,:,:);
   temp8=IRF{kk,1}(Nregion+8,8,:,:);
   temp9=IRF{kk,1}(Nregion+9,9,:,:);
   temp10=IRF{kk,1}(Nregion+10,10,:,:);
   
   ir1(kk,:,:)=squeeze(temp1);
   ir2(kk,:,:)=squeeze(temp2);
   ir3(kk,:,:)=squeeze(temp3);
   ir4(kk,:,:)=squeeze(temp4);
   ir5(kk,:,:)=squeeze(temp5);
   ir6(kk,:,:)=squeeze(temp6);
   ir7(kk,:,:)=squeeze(temp7);
   ir8(kk,:,:)=squeeze(temp8);
   ir9(kk,:,:)=squeeze(temp9);
   ir10(kk,:,:)=squeeze(temp10);
   
   fe1=zeros(1,Nregion,3);
   fe2=zeros(1,Nregion,3);
   fe3=zeros(1,Nregion,3);
   fe4=zeros(1,Nregion,3);
   fe5=zeros(1,Nregion,3);
   fe6=zeros(1,Nregion,3);
   fe7=zeros(1,Nregion,3);
   fe8=zeros(1,Nregion,3);
   fe9=zeros(1,Nregion,3);
   fe10=zeros(1,Nregion,3);
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,1,:,:);    
   fe1(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,2,:,:);    
   fe2(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,3,:,:);    
   fe3(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,4,:,:);    
   fe4(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,5,:,:);    
   fe5(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,6,:,:);    
   fe6(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,7,:,:);    
   fe7(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,8,:,:);    
   fe8(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,9,:,:);    
   fe9(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,10,:,:);    
   fe10(1,jj,:)=squeeze(temp5);    
   end
   
   FE1(kk,:,:)=fe1; FE2(kk,:,:)=fe2; FE3(kk,:,:)=fe3; FE4(kk,:,:)=fe4; 
   FE5(kk,:,:)=fe5; FE6(kk,:,:)=fe6; FE7(kk,:,:)=fe7; FE8(kk,:,:)=fe8; 
   FE9(kk,:,:)=fe9; FE10(kk,:,:)=fe10; 
end

% IRFs of return for region j wrt to illiq shock region j
% do average over 2008 recession.



D1=120:138; 

IR1_1=squeeze(mean(ir1(D1,:,:),1));
IR2_1=squeeze(mean(ir2(D1,:,:),1));
IR3_1=squeeze(mean(ir3(D1,:,:),1));
IR4_1=squeeze(mean(ir4(D1,:,:),1));
IR5_1=squeeze(mean(ir5(D1,:,:),1));
IR6_1=squeeze(mean(ir6(D1,:,:),1));
IR7_1=squeeze(mean(ir7(D1,:,:),1));
IR8_1=squeeze(mean(ir8(D1,:,:),1));
IR9_1=squeeze(mean(ir9(D1,:,:),1));
IR10_1=squeeze(mean(ir10(D1,:,:),1));

figure(1)
subplot(2,5,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
axis tight
subplot(2,5,2)
plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR2_1(:,1),'b--')
plot(0:HO-1,IR2_1(:,3),'b--')
axis tight
subplot(2,5,3)
plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR3_1(:,1),'b--')
plot(0:HO-1,IR3_1(:,3),'b--')
axis tight
subplot(2,5,4)
plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR4_1(:,1),'b--')
plot(0:HO-1,IR4_1(:,3),'b--')
axis tight
subplot(2,5,5)
plot(0:HO-1,IR5_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR5_1(:,1),'b--')
plot(0:HO-1,IR5_1(:,3),'b--')
axis tight
subplot(2,5,6)
plot(0:HO-1,IR6_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR6_1(:,1),'b--')
plot(0:HO-1,IR6_1(:,3),'b--')
axis tight
subplot(2,5,7)
plot(0:HO-1,IR7_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR7_1(:,1),'b--')
plot(0:HO-1,IR7_1(:,3),'b--')
axis tight
subplot(2,5,8)
plot(0:HO-1,IR8_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR8_1(:,1),'b--')
plot(0:HO-1,IR8_1(:,3),'b--')
axis tight
subplot(2,5,9)
plot(0:HO-1,IR9_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR9_1(:,1),'b--')
plot(0:HO-1,IR9_1(:,3),'b--')
axis tight
subplot(2,5,10)
plot(0:HO-1,IR10_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR10_1(:,1),'b--')
plot(0:HO-1,IR10_1(:,3),'b--')
axis tight
matlab2tikz('UK_IRF_RTOV_20.tex')

figure(2)
set(0,'DefaultAxesColorOrder',autumn(10));
subplot(5,2,1)
area(yearlab, 100*FE1(:,:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('EE')
subplot(5,2,2)
area(yearlab, 100*FE2(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
title('EM')
subplot(5,2,3)
area(yearlab, 100*FE3(:,:,2))
axis tight
ylim([0 50])
shadenber()
ylabel('%')
title('NE')
subplot(5,2,4)
area(yearlab, 100*FE4(:,:,2))
axis tight
ylim([0 50])
shadenber()
title('NW')
subplot(5,2,5)
area(yearlab, 100*FE5(:,:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('LO')
subplot(5,2,6)
area(yearlab, 100*FE6(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
title('SE')
subplot(5,2,7)
area(yearlab, 100*FE7(:,:,2))
axis tight
ylim([0 50])
shadenber()
ylabel('%')
title('SW')
subplot(5,2,8)
area(yearlab, 100*FE8(:,:,2))
axis tight
ylim([0 50])
shadenber()
title('WA')
subplot(5,2,9)
area(yearlab, 100*FE7(:,:,2))
axis tight
ylim([0 50])
shadenber()
ylabel('%')
title('West Mid')
subplot(5,2,10)
area(yearlab, 100*FE8(:,:,2))
axis tight
ylim([0 50])
shadenber()
title('YO')
xlabel('Time')
legend('EE','EM','NE','NW','LO','SE','SW','WA','WM','YO','Location','SouthOutside')
matlab2tikz('UK_FEVD_ILLIQ_RTOV_R.tex')
%
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
matlab2tikz('UK_TIC_RtoV_R.tex')

% net directional regions returns
figure(12)
subplot(2,5,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('EE')
subplot(2,5,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('EM')
subplot(2,5,3)
plot(yearlab,NDC(:,3,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,3,1),'k--')
plot(yearlab,NDC(:,3,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('NE')
subplot(2,5,4)
plot(yearlab,NDC(:,4,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,4,1),'k--')
plot(yearlab,NDC(:,4,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('NW')
subplot(2,5,5)
plot(yearlab,NDC(:,5,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,5,1),'k--')
plot(yearlab,NDC(:,5,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('LO')
subplot(2,5,6)
plot(yearlab,NDC(:,6,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,6,1),'k--')
plot(yearlab,NDC(:,6,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('SE')
subplot(2,5,7)
plot(yearlab,NDC(:,7,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,7,1),'k--')
plot(yearlab,NDC(:,7,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('SW')
subplot(2,5,8)
plot(yearlab,NDC(:,8,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,8,1),'k--')
plot(yearlab,NDC(:,8,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('WA')
subplot(2,5,9)
plot(yearlab,NDC(:,9,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,9,1),'k--')
plot(yearlab,NDC(:,9,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('WM')
subplot(2,5,10)
plot(yearlab,NDC(:,10,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,10,1),'k--')
plot(yearlab,NDC(:,10,3),'k--')
axis tight
ylim([-1 3])
shadenber()
title('YO')
matlab2tikz('UK_NDC_RTV_RET_R.tex')


% net directional regions illiq
figure(13)
subplot(2,5,1)
plot(yearlab,NDC(:,11,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,11,1),'k--')
plot(yearlab,NDC(:,11,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('EE')
subplot(2,5,2)
plot(yearlab,NDC(:,12,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,12,1),'k--')
plot(yearlab,NDC(:,12,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('EM')
subplot(2,5,3)
plot(yearlab,NDC(:,13,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,13,1),'k--')
plot(yearlab,NDC(:,13,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('NE')
subplot(2,5,4)
plot(yearlab,NDC(:,14,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,14,1),'k--')
plot(yearlab,NDC(:,14,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('NW')
subplot(2,5,5)
plot(yearlab,NDC(:,15,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,15,1),'k--')
plot(yearlab,NDC(:,15,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('LO')
subplot(2,5,6)
plot(yearlab,NDC(:,16,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,16,1),'k--')
plot(yearlab,NDC(:,16,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('SE')
subplot(2,5,7)
plot(yearlab,NDC(:,17,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,17,1),'k--')
plot(yearlab,NDC(:,17,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('SW')
subplot(2,5,8)
plot(yearlab,NDC(:,18,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,18,1),'k--')
plot(yearlab,NDC(:,18,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('WA')
subplot(2,5,9)
plot(yearlab,NDC(:,19,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,19,1),'k--')
plot(yearlab,NDC(:,19,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('WM')
subplot(2,5,10)
plot(yearlab,NDC(:,20,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,20,1),'k--')
plot(yearlab,NDC(:,20,3),'k--')
axis tight
ylim([-1.5 2.5])
shadenber()
title('YO')
matlab2tikz('UK_NDC_RTV_ILQ_R.tex')



%%

% REGIONAL IVOL NOW

load R_IV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);
ir2=ir1; ir3=ir1; ir4=ir1; ir5=ir1; ir6=ir1; ir7=ir1; ir8=ir1; ir9=ir1; ir10=ir1;

FE1=zeros(TTT,Nregion,3);
FE2=FE1; FE3=FE1; FE4=FE1; FE5=FE1; FE6=FE1; FE7=FE1; FE8=FE1; FE9=FE1; FE10=FE1;

parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);
   temp2=IRF{kk,1}(Nregion+2,2,:,:);
   temp3=IRF{kk,1}(Nregion+3,3,:,:);
   temp4=IRF{kk,1}(Nregion+4,4,:,:);
   temp5=IRF{kk,1}(Nregion+1,5,:,:);
   temp6=IRF{kk,1}(Nregion+2,6,:,:);
   temp7=IRF{kk,1}(Nregion+3,7,:,:);
   temp8=IRF{kk,1}(Nregion+4,8,:,:);
   temp9=IRF{kk,1}(Nregion+3,9,:,:);
   temp10=IRF{kk,1}(Nregion+4,10,:,:);
   
   ir1(kk,:,:)=squeeze(temp1);
   ir2(kk,:,:)=squeeze(temp2);
   ir3(kk,:,:)=squeeze(temp3);
   ir4(kk,:,:)=squeeze(temp4);
   ir5(kk,:,:)=squeeze(temp5);
   ir6(kk,:,:)=squeeze(temp6);
   ir7(kk,:,:)=squeeze(temp7);
   ir8(kk,:,:)=squeeze(temp8);
   ir9(kk,:,:)=squeeze(temp9);
   ir10(kk,:,:)=squeeze(temp10);
   
   fe1=zeros(1,Nregion,3);
   fe2=zeros(1,Nregion,3);
   fe3=zeros(1,Nregion,3);
   fe4=zeros(1,Nregion,3);
   fe5=zeros(1,Nregion,3);
   fe6=zeros(1,Nregion,3);
   fe7=zeros(1,Nregion,3);
   fe8=zeros(1,Nregion,3);
   fe9=zeros(1,Nregion,3);
   fe10=zeros(1,Nregion,3);
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,1,:,:);    
   fe1(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,2,:,:);    
   fe2(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,3,:,:);    
   fe3(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,4,:,:);    
   fe4(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,5,:,:);    
   fe5(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,6,:,:);    
   fe6(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,7,:,:);    
   fe7(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,8,:,:);    
   fe8(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,9,:,:);    
   fe9(1,jj,:)=squeeze(temp5);    
   end
   
   for jj=1:Nregion
   temp5=FEV{kk,1}(Nregion+jj,10,:,:);    
   fe10(1,jj,:)=squeeze(temp5);    
   end
   
   FE1(kk,:,:)=fe1; FE2(kk,:,:)=fe2; FE3(kk,:,:)=fe3; FE4(kk,:,:)=fe4; 
   FE5(kk,:,:)=fe5; FE6(kk,:,:)=fe6; FE7(kk,:,:)=fe7; FE8(kk,:,:)=fe8; 
   FE9(kk,:,:)=fe9; FE10(kk,:,:)=fe10; 
end

% IRFs of return for region j wrt to illiq shock region j
% do average over 2008 recession.



D1=120:138; 

IR1_1=squeeze(mean(ir1(D1,:,:),1));
IR2_1=squeeze(mean(ir2(D1,:,:),1));
IR3_1=squeeze(mean(ir3(D1,:,:),1));
IR4_1=squeeze(mean(ir4(D1,:,:),1));
IR5_1=squeeze(mean(ir5(D1,:,:),1));
IR6_1=squeeze(mean(ir6(D1,:,:),1));
IR7_1=squeeze(mean(ir7(D1,:,:),1));
IR8_1=squeeze(mean(ir8(D1,:,:),1));
IR9_1=squeeze(mean(ir9(D1,:,:),1));
IR10_1=squeeze(mean(ir10(D1,:,:),1));

figure(10)
subplot(2,5,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
axis tight
title('EE')
subplot(2,5,2)
plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR2_1(:,1),'b--')
plot(0:HO-1,IR2_1(:,3),'b--')
axis tight
title('EM')
subplot(2,5,3)
plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR3_1(:,1),'b--')
plot(0:HO-1,IR3_1(:,3),'b--')
axis tight
title('NE')
subplot(2,5,4)
plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR4_1(:,1),'b--')
plot(0:HO-1,IR4_1(:,3),'b--')
axis tight
title('NW')
subplot(2,5,5)
plot(0:HO-1,IR5_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR5_1(:,1),'b--')
plot(0:HO-1,IR5_1(:,3),'b--')
axis tight
title('LO')
subplot(2,5,6)
plot(0:HO-1,IR6_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR6_1(:,1),'b--')
plot(0:HO-1,IR6_1(:,3),'b--')
axis tight
title('SE')
subplot(2,5,7)
plot(0:HO-1,IR7_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR7_1(:,1),'b--')
plot(0:HO-1,IR7_1(:,3),'b--')
axis tight
title('SW')
subplot(2,5,8)
plot(0:HO-1,IR8_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR8_1(:,1),'b--')
plot(0:HO-1,IR8_1(:,3),'b--')
axis tight
title('WA')
subplot(2,5,9)
plot(0:HO-1,IR9_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR9_1(:,1),'b--')
plot(0:HO-1,IR9_1(:,3),'b--')
axis tight
title('WM')
subplot(2,5,10)
plot(0:HO-1,IR10_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR10_1(:,1),'b--')
plot(0:HO-1,IR10_1(:,3),'b--')
axis tight
title('YO')
matlab2tikz('UK_IRF_invV_20.tex')

figure(999)
subplot(5,2,1)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
title('EE')
subplot(5,2,2)
mesh(yearlab,0:HO-1,ir2(:,:,2)')
colormap autumn
axis tight
title('EM')
subplot(5,2,3)
mesh(yearlab,0:HO-1,ir3(:,:,2)')
colormap autumn
axis tight
title('LO')
subplot(5,2,4)
mesh(yearlab,0:HO-1,ir4(:,:,2)')
colormap autumn
axis tight
title('NE')
subplot(5,2,5)
mesh(yearlab,0:HO-1,ir5(:,:,2)')
colormap autumn
axis tight
title('NW')
subplot(5,2,6)
mesh(yearlab,0:HO-1,ir6(:,:,2)')
colormap autumn
axis tight
title('SE')
subplot(5,2,7)
mesh(yearlab,0:HO-1,ir7(:,:,2)')
colormap autumn
axis tight
title('SW')
subplot(5,2,8)
mesh(yearlab,0:HO-1,ir8(:,:,2)')
colormap autumn
axis tight
title('WA')
subplot(5,2,9)
mesh(yearlab,0:HO-1,ir9(:,:,2)')
colormap autumn
axis tight
title('SW')
subplot(5,2,10)
mesh(yearlab,0:HO-1,ir10(:,:,2)')
colormap autumn
axis tight
title('YO')
matlab2tikz('UK_3D_R.tex')


figure(20)
set(0,'DefaultAxesColorOrder',autumn(10));
subplot(5,2,1)
area(yearlab, 100*FE1(:,:,2))
axis tight
ylabel('%')
ylim([0 30])
shadenber()
title('EE')
subplot(5,2,2)
area(yearlab, 100*FE2(:,:,2))
axis tight
ylim([0 30])
shadenber()
xlabel('Time')
title('EM')
subplot(5,2,3)
area(yearlab, 100*FE3(:,:,2))
axis tight
ylim([0 30])
shadenber()
ylabel('%')
title('NE')
subplot(5,2,4)
area(yearlab, 100*FE4(:,:,2))
axis tight
ylim([0 30])
shadenber()
title('NW')
subplot(5,2,5)
area(yearlab, 100*FE5(:,:,2))
axis tight
ylabel('%')
ylim([0 30])
shadenber()
title('LO')
subplot(5,2,6)
area(yearlab, 100*FE6(:,:,2))
axis tight
ylim([0 30])
shadenber()
xlabel('Time')
title('SE')
subplot(5,2,7)
area(yearlab, 100*FE7(:,:,2))
axis tight
ylim([0 30])
shadenber()
ylabel('%')
title('SW')
subplot(5,2,8)
area(yearlab, 100*FE8(:,:,2))
axis tight
ylim([0 30])
shadenber()
title('WA')
subplot(5,2,9)
area(yearlab, 100*FE7(:,:,2))
axis tight
ylim([0 30])
shadenber()
ylabel('%')
title('West Mid')
subplot(5,2,10)
area(yearlab, 100*FE8(:,:,2))
axis tight
ylim([0 30])
shadenber()
title('YO')
xlabel('Time')
legend('EE','EM','NE','NW','LO','SE','SW','WA','WM','YO','Location','SouthOutside')
matlab2tikz('UK_FEVD_ILLIQ_invV_R.tex')


%%
% Now for network stuff

% total connectedness
figure(110)
plot(yearlab,TIC(:,2),'b-','LineWidth',1.1)
hold on,
plot(yearlab,TIC(:,1),'b--')
plot(yearlab,TIC(:,3),'b--')
shadenber()
axis tight
xlabel('Time')
matlab2tikz('UK_TIC_invV_R.tex')

% net directional regions returns
figure(120)
subplot(2,5,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('EE')
subplot(2,5,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('EM')
subplot(2,5,3)
plot(yearlab,NDC(:,3,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,3,1),'k--')
plot(yearlab,NDC(:,3,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('NE')
subplot(2,5,4)
plot(yearlab,NDC(:,4,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,4,1),'k--')
plot(yearlab,NDC(:,4,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('NW')
subplot(2,5,5)
plot(yearlab,NDC(:,5,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,5,1),'k--')
plot(yearlab,NDC(:,5,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('LO')
subplot(2,5,6)
plot(yearlab,NDC(:,6,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,6,1),'k--')
plot(yearlab,NDC(:,6,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('SE')
subplot(2,5,7)
plot(yearlab,NDC(:,7,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,7,1),'k--')
plot(yearlab,NDC(:,7,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('SW')
subplot(2,5,8)
plot(yearlab,NDC(:,8,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,8,1),'k--')
plot(yearlab,NDC(:,8,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('WA')
subplot(2,5,9)
plot(yearlab,NDC(:,9,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,9,1),'k--')
plot(yearlab,NDC(:,9,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('WM')
subplot(2,5,10)
plot(yearlab,NDC(:,10,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,10,1),'k--')
plot(yearlab,NDC(:,10,3),'k--')
axis tight
ylim([-2.5 1.5])
shadenber()
title('YO')
matlab2tikz('UK_NDC_IV_RET_R.tex')


% net directional regions illiq
figure(130)
subplot(2,5,1)
plot(yearlab,NDC(:,11,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,11,1),'k--')
plot(yearlab,NDC(:,11,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('EE')
subplot(2,5,2)
plot(yearlab,NDC(:,12,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,12,1),'k--')
plot(yearlab,NDC(:,12,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('EM')
subplot(2,5,3)
plot(yearlab,NDC(:,13,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,13,1),'k--')
plot(yearlab,NDC(:,13,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('NE')
subplot(2,5,4)
plot(yearlab,NDC(:,14,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,14,1),'k--')
plot(yearlab,NDC(:,14,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('NW')
subplot(2,5,5)
plot(yearlab,NDC(:,15,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,15,1),'k--')
plot(yearlab,NDC(:,15,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('LO')
subplot(2,5,6)
plot(yearlab,NDC(:,16,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,16,1),'k--')
plot(yearlab,NDC(:,16,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('SE')
subplot(2,5,7)
plot(yearlab,NDC(:,17,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,17,1),'k--')
plot(yearlab,NDC(:,17,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('SW')
subplot(2,5,8)
plot(yearlab,NDC(:,18,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,18,1),'k--')
plot(yearlab,NDC(:,18,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('WA')
subplot(2,5,9)
plot(yearlab,NDC(:,19,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,19,1),'k--')
plot(yearlab,NDC(:,19,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('WM')
subplot(2,5,10)
plot(yearlab,NDC(:,20,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,20,1),'k--')
plot(yearlab,NDC(:,20,3),'k--')
axis tight
ylim([-0.5 3.5])
shadenber()
title('YO')
matlab2tikz('UK_NDC_IV_ILQ_R.tex')