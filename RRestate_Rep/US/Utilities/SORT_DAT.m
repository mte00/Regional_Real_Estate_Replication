% SORTS OUT DATA
LIQM=data(:,3:5);
% RESCALE 1/V and RtoV
LIQM(:,1)=LIQM(:,1)*10^8; LIQM(:,3)=LIQM(:,3)*10^2;
% REARRANGE PROXY TO CORRESPOND TO 1=ROLL, 2=1/V, 3=RtoV.
LIQM=[LIQM(:,2), LIQM(:,1), LIQM(:,3)];
% CREATE TS MOVING AVERAGE
LIQMA=tsmovavg(LIQM(:,2:end),'s',4,1);
%LIQMA=movavg(LIQM(:,2:end),'simple',4,'fill');    
LIQM2=100*(LIQM(:,2:end)-LIQMA)./LIQMA;    
LIQM=[LIQM(:,1),LIQM2];
% GET OTHER DATA and sort into following order: house price, u, w, s, c, p
data=[data(:,2), data(:,1), data(:,6:end)];
% GENERATE TRANSFORMATIONS
% Quarterly growth in data apart from last column which is annual growth
% models span from 1998Q1--2017Q3 79 observations in total
for i=1:size(data,2)-1
    if i==2
        temp=(1+data(2:end,i)./100).^(1/4)-1;
        temp=100*temp;
    else
        temp=100*log(data(2:end,i)./data(1:end-1,i));
    end
    data2(:,i)=temp;
end
data3=100*log(data(5:end,end)./data(1:end-4,end));
data2=data2(4:end,:);
% NOW LOAD IN MORTGAGE RATE AND CREDIT DATA
rr=xlsread('mortgage_rate','Sheet1','B13:B91');
C=xlsread('credit','Sheet1','E21:E99');

data=[LIQM(5:end,ILLIQ), data2,data3,rr,C];

