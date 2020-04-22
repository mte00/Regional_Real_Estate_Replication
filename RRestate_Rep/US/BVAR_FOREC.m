function [out] = BVAR_FOREC(data,L,nsim,shrinkage,T_thres,nfore,fdirname)

[T,N]=size(data); K=N*L+1;
% Get data in VAR setup
dat2=data;
for i=1:L
    temp=lag0(dat2,i);
    X(:,1+N*(i-1):i*N)=temp(1+L:T,:);
end
y=dat2(1+L:T,:); T=T-L; X=[ones(T,1),X];
clear temp
T_full=T;
anum=T_full-T_thres+1;
PL_ave=single(zeros(anum,nfore));
y_foredraws=single(zeros(nsim,nfore,N));
PL_draws=single(zeros(nsim,nfore));
%%
file=0;
tic
for sample=T_thres:T_full
    disp(['Now you are running sample ' num2str(sample) ' of ' num2str(T_full)] )
    fname=strcat(fdirname,'Forecast',num2str(file),'.mat'); %name of file to save forecast density
   if sample<=T_full-nfore
       Y_f=y(sample+1:sample+nfore,:); %OOS observations
   else
       Y_f=[y(sample+1:T_full,:);NaN(nfore-(T_full-sample),N)];
   end
    
   y1=y(1:sample,:);
   X1=X(1:sample,:);
   tt=sample;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SI,PI,a,RI]=Minn_NWprior(y1,tt,N,L,shrinkage);
% SI prior mean matrix. PI prior variance of SI (diagonal matrix (K x K)),
% and RI is prior for covariance matrix of VAR model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate weights
weights=normker(tt,sqrt(tt));
% Follow Petrova (2018) work with precision
priorprec0=PI^(-1);
priorprec0=sparse(priorprec0); % create sparse matrix allowing more efficient allocation of memory
RI=sparse(RI); % create sparse matrix allowing more efficient allocation of memory
% clear priorprec0 
% run computation using standard and then sparse arrays to check for
% consistency and then uncomment clear PI and clear priorprec0
% B=(X'*diag(w)*X)^(-1)*(X'*diag(w)*y);
% bhat=B(:);
% bhat=((kron(eye(N),X'*diag(w)*X))^(-1))*(kron(eye(N),X'*diag(w))*y(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posteriors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   bayesprec=(priorprec0+X1'*X1);
   bayessv=bayesprec^(-1);
   BB=bayessv*(X1'*y1+priorprec0*SI);
   bayesb=BB(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  bayesalpha=a+tt;
  g1=SI'*priorprec0*SI;
  g2=y1'*y1;
  g3=BB'*bayesprec*BB;
  bayesgamma=RI+g1+g2-g3;
  bayesgamma=0.5*bayesgamma+0.5*bayesgamma'; %it is symmetric but just in case
  %make draws
  parfor ii=1:nsim
     mm=0;
     while mm<1
         SIGMA=single(iwishrnd(bayesgamma,bayesalpha)); % Draw from IW distribution
         % Fi=mvnrnd(bayesb,kron(Sigma,bayesv))';
         nu=randn(N*L+1,N);
         Fi1=single((BB+chol(bayessv)'*nu*(chol(SIGMA)))');
%         Ficom(1:N,:)=(Fi1(:,2:N*L+1));
%         for pp=2:L
%            Ficom(1+N*(pp-1):pp*N,1+N*(pp-2):N*(pp-1))=eye(N); % companion without constant 
%         end
         max_eig=max(abs(eig([Fi1(:,2:end);eye(N),zeros(N,N)])));
         if max_eig<.999 % check stability of draw
             stab_ind=1;
             mm=1;
         end
     end
     
%     if kk==tt
        % Do Forecasting
        ylast=y1(tt,:);
        xlast=X1(tt,2:N*(L-1)+1);
%        miu=[Fi(:,1); zeros(N*(L-1),1)];
%        BB1=Fi(:,2:end);
%        BETA_FORE=[BB1; eye(N*(L-1)) zeros(N*(L-1),N)];
        y_fore=zeros(nfore,N);
%        VM=0;
        for pp=1:nfore
%           VM=VM+(BETA_FORE^(pp-1))*miu;
%           FORECASTS=VM+(BETA_FORE^(pp))*X_f';
%           y_fore1(pp,:)=FORECASTS(1:N,:)'; % THIS IS ALL FOR DIRECT
%           FORECASTING AND WE ARE USING ITERATIVE.
            xlast=[ylast, xlast(1:N*(L-1))];
            uu=randn(N,1);
            ylast=(Fi1*[1, xlast]'+(uu'*chol(SIGMA))')';
            y_fore(pp,:)=ylast;
           PL_draws(ii,pp)=mvnpdf(Y_f(pp,:),y_fore(pp,:),SIGMA);
        end
        y_foredraws(ii,:,:)=y_fore;
%     end
   end


Y_fore=y_foredraws; PL_full=PL_draws;

if sample==T_full
    
for jj=1:nfore
    PL_ave(sample-T_thres+1,jj)=mean(PL_draws(:,jj));
end

save(fname,'Y_fore','PL_full','PL_ave');

else

save(fname,'Y_fore','PL_full');
end


file=file+1;
end
out=1;
toc

end