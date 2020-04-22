clear; clc;
addpath('src')
addpath('Results')
addpath('Utilities')

HO=21; % horizon we want to look at for FEVDS
Nregion=4; % 4 regions for US in our study
% first do RtoV regional
yearlab=(1985.00:(1/12):2018.75)';

load R_RtoV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);
ir2=ir1; ir3=ir1; ir4=ir1;

FE1=zeros(TTT,Nregion,3);
FE2=FE1; FE3=FE1; FE4=FE1;

parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);
   temp2=IRF{kk,1}(Nregion+2,2,:,:);
   temp3=IRF{kk,1}(Nregion+3,3,:,:);
   temp4=IRF{kk,1}(Nregion+4,4,:,:);
   ir1(kk,:,:)=squeeze(temp1);
   ir2(kk,:,:)=squeeze(temp2);
   ir3(kk,:,:)=squeeze(temp3);
   ir4(kk,:,:)=squeeze(temp4);
   
   fe1=zeros(1,Nregion,3);
   fe2=zeros(1,Nregion,3);
   fe3=zeros(1,Nregion,3);
   fe4=zeros(1,Nregion,3);
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
   
   FE1(kk,:,:)=fe1; FE2(kk,:,:)=fe2; FE3(kk,:,:)=fe3; FE4(kk,:,:)=fe4; 
end

% IRFs of return for region j wrt to illiq shock region j

%D1=1:239; D2=240:406;
%IR1_1=squeeze(mean(ir1(D1,:,:),1));
%IR1_2=squeeze(mean(ir1(D2,:,:),1));
%IR2_1=squeeze(mean(ir2(D1,:,:),1));
%IR2_2=squeeze(mean(ir2(D2,:,:),1));
%IR3_1=squeeze(mean(ir3(D1,:,:),1));
%IR3_2=squeeze(mean(ir3(D2,:,:),1));
%IR4_1=squeeze(mean(ir4(D1,:,:),1));
%IR4_2=squeeze(mean(ir4(D2,:,:),1));


%figure(1)
%subplot(2,4,1)
%plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR1_1(:,1),'b--')
%plot(0:HO-1,IR1_1(:,3),'b--')
%title('Midwest')
%ylabel('% 1985-2004')
%subplot(2,4,2)
%plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR2_1(:,1),'b--')
%plot(0:HO-1,IR2_1(:,3),'b--')
%title('Northeast')
%subplot(2,4,3)
%plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR3_1(:,1),'b--')
%plot(0:HO-1,IR3_1(:,3),'b--')
%title('South')
%subplot(2,4,4)
%plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR4_1(:,1),'b--')
%plot(0:HO-1,IR4_1(:,3),'b--')
%title('West')
%subplot(2,4,5)
%plot(0:HO-1,IR1_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR1_2(:,1),'b--')
%plot(0:HO-1,IR1_2(:,3),'b--')
%ylabel('% 2005-2018')
%xlabel('Horizon')
%subplot(2,4,6)
%plot(0:HO-1,IR2_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR2_2(:,1),'b--')
%plot(0:HO-1,IR2_2(:,3),'b--')
%xlabel('Horizon')
%subplot(2,4,7)
%plot(0:HO-1,IR3_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR3_2(:,1),'b--')
%plot(0:HO-1,IR3_2(:,3),'b--')
%xlabel('Horizon')
%subplot(2,4,8)
%plot(0:HO-1,IR4_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR4_2(:,1),'b--')
%plot(0:HO-1,IR4_2(:,3),'b--')
%xlabel('Horizon')
%matlab2tikz('IRF_RTOV_20.tex')
D1=276:294;
IR1_1=squeeze(mean(ir1(D1,:,:),1));
IR2_1=squeeze(mean(ir2(D1,:,:),1));
IR3_1=squeeze(mean(ir3(D1,:,:),1));
IR4_1=squeeze(mean(ir4(D1,:,:),1));
figure(1)
subplot(2,2,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
axis tight
title('MW')
subplot(2,2,2)
plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR2_1(:,1),'b--')
plot(0:HO-1,IR2_1(:,3),'b--')
axis tight
title('NE')
subplot(2,2,3)
plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR3_1(:,1),'b--')
plot(0:HO-1,IR3_1(:,3),'b--')
axis tight
title('S')
subplot(2,2,4)
plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR4_1(:,1),'b--')
plot(0:HO-1,IR4_1(:,3),'b--')
axis tight
title('W')
matlab2tikz('US_IRF_RTOV_20.tex')

figure(999)
subplot(2,2,1)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
title('MW')
subplot(2,2,2)
mesh(yearlab,0:HO-1,ir2(:,:,2)')
colormap autumn
axis tight
title('NE')
subplot(2,2,3)
mesh(yearlab,0:HO-1,ir3(:,:,2)')
colormap autumn
axis tight
title('W')
subplot(2,2,4)
mesh(yearlab,0:HO-1,ir4(:,:,2)')
colormap autumn
axis tight
title('S')
matlab2tikz('US_3D_R.tex')

figure(2)
set(0,'DefaultAxesColorOrder',autumn(4));
subplot(2,2,1)
area(yearlab, 100*FE1(:,:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('MW')
subplot(2,2,2)
area(yearlab, 100*FE2(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
title('NE')
subplot(2,2,3)
area(yearlab, 100*FE3(:,:,2))
axis tight
ylim([0 50])
shadenber()
ylabel('%')
title('S')
subplot(2,2,4)
area(yearlab, 100*FE4(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
legend('MW','NE','S','W','Location','SouthOutside')
title('W')
matlab2tikz('US_FEVD_ILLIQ_RTOV_R.tex')
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
matlab2tikz('US_TIC_RtoV_R.tex')

% net directional regions returns
figure(12)
subplot(2,2,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('MW')
subplot(2,2,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('NE')
subplot(2,2,3)
plot(yearlab,NDC(:,3,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,3,1),'k--')
plot(yearlab,NDC(:,3,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('S')
subplot(2,2,4)
plot(yearlab,NDC(:,4,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,4,1),'k--')
plot(yearlab,NDC(:,4,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('W')
matlab2tikz('US_NDC_RTV_RET_R.tex')


% net directional regions illiq
figure(13)
subplot(2,2,1)
plot(yearlab,NDC(:,5,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,5,1),'k--')
plot(yearlab,NDC(:,5,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('MW')
subplot(2,2,2)
plot(yearlab,NDC(:,6,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,6,1),'k--')
plot(yearlab,NDC(:,6,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('NE')
subplot(2,2,3)
plot(yearlab,NDC(:,7,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,7,1),'k--')
plot(yearlab,NDC(:,7,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('S')
subplot(2,2,4)
plot(yearlab,NDC(:,8,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,8,1),'k--')
plot(yearlab,NDC(:,8,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('W')
matlab2tikz('US_NDC_RTV_ILQ_R.tex')



%%

% REGIONAL IVOL NOW

load R_IV_gr.mat

% get impulse response functions of return wrt its corresponding illiq
% measure.

ir1=zeros(TTT,HO,3);
ir2=ir1; ir3=ir1; ir4=ir1;

FE1=zeros(TTT,Nregion,3);
FE2=FE1; FE3=FE1; FE4=FE1;

parfor kk=1:TTT
   temp1=IRF{kk,1}(Nregion+1,1,:,:);
   temp2=IRF{kk,1}(Nregion+2,2,:,:);
   temp3=IRF{kk,1}(Nregion+3,3,:,:);
   temp4=IRF{kk,1}(Nregion+4,4,:,:);
   ir1(kk,:,:)=squeeze(temp1);
   ir2(kk,:,:)=squeeze(temp2);
   ir3(kk,:,:)=squeeze(temp3);
   ir4(kk,:,:)=squeeze(temp4);
   
   fe1=zeros(1,Nregion,3);
   fe2=zeros(1,Nregion,3);
   fe3=zeros(1,Nregion,3);
   fe4=zeros(1,Nregion,3);
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
   
   FE1(kk,:,:)=fe1; FE2(kk,:,:)=fe2; FE3(kk,:,:)=fe3; FE4(kk,:,:)=fe4; 
end

% IRFs of return for region j wrt to illiq shock region j

%D1=1:239; D2=240:406;
%IR1_1=squeeze(mean(ir1(D1,:,:),1));
%IR1_2=squeeze(mean(ir1(D2,:,:),1));
%IR2_1=squeeze(mean(ir2(D1,:,:),1));
%IR2_2=squeeze(mean(ir2(D2,:,:),1));
%IR3_1=squeeze(mean(ir3(D1,:,:),1));
%IR3_2=squeeze(mean(ir3(D2,:,:),1));
%IR4_1=squeeze(mean(ir4(D1,:,:),1));
%IR4_2=squeeze(mean(ir4(D2,:,:),1));


%figure(1)
%subplot(2,4,1)
%plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR1_1(:,1),'b--')
%plot(0:HO-1,IR1_1(:,3),'b--')
%title('Midwest')
%ylabel('% 1985-2004')
%subplot(2,4,2)
%plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR2_1(:,1),'b--')
%plot(0:HO-1,IR2_1(:,3),'b--')
%title('Northeast')
%subplot(2,4,3)
%plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR3_1(:,1),'b--')
%plot(0:HO-1,IR3_1(:,3),'b--')
%title('South')
%subplot(2,4,4)
%plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR4_1(:,1),'b--')
%plot(0:HO-1,IR4_1(:,3),'b--')
%title('West')
%subplot(2,4,5)
%plot(0:HO-1,IR1_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR1_2(:,1),'b--')
%plot(0:HO-1,IR1_2(:,3),'b--')
%ylabel('% 2005-2018')
%xlabel('Horizon')
%subplot(2,4,6)
%plot(0:HO-1,IR2_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR2_2(:,1),'b--')
%plot(0:HO-1,IR2_2(:,3),'b--')
%xlabel('Horizon')
%subplot(2,4,7)
%plot(0:HO-1,IR3_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR3_2(:,1),'b--')
%plot(0:HO-1,IR3_2(:,3),'b--')
%xlabel('Horizon')
%subplot(2,4,8)
%plot(0:HO-1,IR4_2(:,2),'b-','LineWidth',1.1)
%hold on,
%plot(0:HO-1,IR4_2(:,1),'b--')
%plot(0:HO-1,IR4_2(:,3),'b--')
%xlabel('Horizon')
%matlab2tikz('IRF_IV_20.tex')
D1=276:294;
IR1_1=squeeze(mean(ir1(D1,:,:),1));
IR2_1=squeeze(mean(ir2(D1,:,:),1));
IR3_1=squeeze(mean(ir3(D1,:,:),1));
IR4_1=squeeze(mean(ir4(D1,:,:),1));
figure(3)
subplot(2,2,1)
plot(0:HO-1,IR1_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR1_1(:,1),'b--')
plot(0:HO-1,IR1_1(:,3),'b--')
axis tight
title('MW')
subplot(2,2,2)
plot(0:HO-1,IR2_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR2_1(:,1),'b--')
plot(0:HO-1,IR2_1(:,3),'b--')
axis tight
title('NE')
subplot(2,2,3)
plot(0:HO-1,IR3_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR3_1(:,1),'b--')
plot(0:HO-1,IR3_1(:,3),'b--')
axis tight
title('S')
subplot(2,2,4)
plot(0:HO-1,IR4_1(:,2),'b-','LineWidth',1.1)
hold on,
plot(0:HO-1,IR4_1(:,1),'b--')
plot(0:HO-1,IR4_1(:,3),'b--')
axis tight
title('W')
matlab2tikz('US_IRF_IV_20.tex')

figure(999)
subplot(2,2,1)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
title('MW')
subplot(2,2,2)
mesh(yearlab,0:HO-1,ir2(:,:,2)')
colormap autumn
axis tight
title('NE')
subplot(2,2,3)
mesh(yearlab,0:HO-1,ir3(:,:,2)')
colormap autumn
axis tight
title('W')
subplot(2,2,4)
mesh(yearlab,0:HO-1,ir1(:,:,2)')
colormap autumn
axis tight
title('S')
matlab2tikz('US_3D_R.tex')

figure(4)
set(0,'DefaultAxesColorOrder',autumn(4));
subplot(2,2,1)
area(yearlab, 100*FE1(:,:,2))
axis tight
ylabel('%')
ylim([0 50])
shadenber()
title('MW')
subplot(2,2,2)
area(yearlab, 100*FE2(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
title('NE')
subplot(2,2,3)
area(yearlab, 100*FE3(:,:,2))
axis tight
ylim([0 50])
shadenber()
ylabel('%')
title('S')
subplot(2,2,4)
area(yearlab, 100*FE4(:,:,2))
axis tight
ylim([0 50])
shadenber()
xlabel('Time')
legend('MW','NE','S','W','Location','SouthOutside')
title('W')
matlab2tikz('US_FEVD_ILLIQ_IV_R.tex')

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
matlab2tikz('US_TIC_RIV_R.tex')

% net directional regions returns
figure(22)
subplot(2,2,1)
plot(yearlab,NDC(:,1,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,1,1),'k--')
plot(yearlab,NDC(:,1,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('MW')
subplot(2,2,2)
plot(yearlab,NDC(:,2,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,2,1),'k--')
plot(yearlab,NDC(:,2,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('NE')
subplot(2,2,3)
plot(yearlab,NDC(:,3,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,3,1),'k--')
plot(yearlab,NDC(:,3,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('S')
subplot(2,2,4)
plot(yearlab,NDC(:,4,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,4,1),'k--')
plot(yearlab,NDC(:,4,3),'k--')
axis tight
ylim([-3 8])
shadenber()
title('W')
matlab2tikz('US_NDC_IV_RET_R.tex')

% net directional regions illiq
figure(23)
subplot(2,2,1)
plot(yearlab,NDC(:,5,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,5,1),'k--')
plot(yearlab,NDC(:,5,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('MW')
subplot(2,2,2)
plot(yearlab,NDC(:,6,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,6,1),'k--')
plot(yearlab,NDC(:,6,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('NE')
subplot(2,2,3)
plot(yearlab,NDC(:,7,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,7,1),'k--')
plot(yearlab,NDC(:,7,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('S')
subplot(2,2,4)
plot(yearlab,NDC(:,8,2),'k-','LineWidth',1.1)
hold on,
plot(yearlab,NDC(:,8,1),'k--')
plot(yearlab,NDC(:,8,3),'k--')
axis tight
ylim([-2 2])
shadenber()
title('W')
matlab2tikz('US_NDC_IV_ILQ_R.tex')
