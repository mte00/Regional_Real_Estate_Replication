[T,N]=size(data);
nsim=5000;
shrinkage=0.05; % Overall shrinkage of Minnesota Prior. Change this to investigate influence
L=2; % 2 lags corresponds to 2 months
K=N*L+1;
% Get data in VAR setup
dat2=data;
for i=1:L
    temp=lag0(dat2,i);
    X(:,1+N*(i-1):i*N)=temp(1+L:T,:);
end
y=dat2(1+L:T,:); T=T-L; X=[ones(T,1),X];
clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SI,PI,a,RI]=Minn_NWprior(dat2,T,N,L,shrinkage);
% SI prior mean matrix. PI prior variance of SI (diagonal matrix (K x K)),
% and RI is prior for covariance matrix of VAR model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate weights
weights=normker(T,sqrt(T));
% Follow Petrova (2018) work with precision
priorprec0=PI^(-1);
clear PI
priorprec0=sparse(priorprec0); % create sparse matrix allowing more efficient allocation of memory
RI=sparse(RI); % create sparse matrix allowing more efficient allocation of memory


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posteriors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stab_ind=zeros(nsim,T,'single');
max_eig=zeros(nsim,T,'single');

% Here i allocate cells to store the time varying coefficient matrices and
% the time varying covariance matrices of the daily data in order to save
% on memory.
% Set variable names:
varname(1,:)='TIC'; varname(2,:)='NDC';varname(3,:)='PWC';
varname(4,:)='IRF'; varname(5,:)='FEV'; varname(6,:)='NNN'; varname(7,:)='TTT';
TTT=T; TT=T; LL=L; NNN=N; NS=nsim;

TIC=single(zeros(TT,3));
NDC=single(zeros(TT,N,3)); % Net-directional Connectedness for each variable
PWC=single(zeros(N,N,3,TT));
HO=20+1;
IRF=cell(TT,1);
FEV=cell(TT,1);

%%
parfor kk=1:T
   w=weights(kk,:);
   bayesprec=(priorprec0+X'*diag(w)*X);
   bayessv=bayesprec^(-1);
   BB=bayessv*((X'*diag(w))*y+priorprec0*SI);
   bayesb=BB(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  bayesalpha=a+sum(w);
  g1=SI'*priorprec0*SI;
  g2=y'*diag(w)*y;
  g3=BB'*bayesprec*BB;
  bayesgamma=RI+g1+g2-g3;
  bayesgamma=0.5*bayesgamma+0.5*bayesgamma'; %it is symmetric but just in case
  %make draws
  FF=single(zeros(N*K,nsim));
  SS=single(zeros(N*N,nsim));
  
  tic=single(zeros(1,NS)); ndc=single(zeros(N,NS)); pw=single(zeros(N,N,NS));
    irfs=single(zeros(N,N,HO,NS)); %storage for irfs
    fevd=single(zeros(N,N,NS)); %storage for fevd
  
  for ii=1:nsim
     mm=0;
     while mm<1
         SIGMA=single(iwishrnd(bayesgamma,bayesalpha)); % Draw from IW distribution
         nu=randn(N*L+1,N);
         Fi1=single((BB+chol(bayessv)'*nu*(chol(SIGMA))))';

         max_eig(ii,kk)=max(abs(eig([Fi1(:,2:end);eye(N),zeros(N,N)])));
         if max_eig(ii,kk)<.999 % check stability of draw
             stab_ind(ii,kk)=1;
             mm=1;
         end
         
       % GET THE STUFF WE CARE ABOUT  
       [irf]=get_GIRF(Fi1,SIGMA,1,LL,HO-1); 
       irfs(:,:,:,ii)=irf;
       [tc,~,~,nd,pwc,fe]=get_timeconnect(N,HO,irf);
       fevd(:,:,ii)=fe; % get FEV @ 20 month horizon
       tic(:,ii)=tc; ndc(:,ii)=nd; pw(:,:,ii)=pwc; 
     end
  end
  irfs=quantile(irfs,[0.16 0.5 0.84],4);
  fevd=quantile(fevd,[0.16 0.5 0.84],3);
  tic=quantile(tic,[0.16 0.5 0.84],2)
  ndc=quantile(ndc,[0.16 0.5 0.84],2);
  pw=quantile(pw,[0.16 0.5 0.84],3);
  TIC(kk,:)=tic;
  NDC(kk,:,:)=ndc;
  PWC(:,:,:,kk)=pw;
  IRF{kk,:}=irfs;
  FEV{kk,:}=fevd;

end
save(DFILE,varname(1,:),varname(2,:),varname(3,:),varname(4,:),varname(5,:),varname(6,:)...
    ,varname(7,:),'-v7.3');
toc
