% ESTIMATE IRFs and FEVDs from LARGE TVP VAR MODEL using QBLL Methods.
N=9;
dimB=N*(1+N*LL);
%%
% GENERATE STORAGE MATRICES AND INCLUDE SOME PRELIMINARIES
HO=21; %IRF Horizon, Also use this for FEVD horizon.

% If problems saving 5D data will split into 4D where we separate matrices
% by number of shocks. When we come to use N=100 this will be a large task.
IRFS=cell(TT,1);
FEVD=cell(TT,1);

% LOOP OVER TIME
kk=1;
while kk<=TT
    irfs=zeros(N,N,HO,NS); %storage for irfs
    fevd=zeros(N,N,HO,NS); %storage for fevd
    THETA=SD{kk,:};
    SIGMA=OM{kk,:};
   for ii = 1:NS
       sd=THETA(:,ii);
       sd=reshape(sd,N,N*LL+1);
       sig=SIGMA(:,ii);
       sig=reshape(sig,N,N);
       % compute irfs
       [irf]=get_GIRF(sd,sig,1,LL,HO-1);
       %[irf]=impulse(sd,chol(sig)',1,LL,HO-1);
       irfs(:,:,:,ii)=irf;
       fev=vardecomp(N,HO,irf);
       fevd(:,:,:,ii)=fev;
   end
   IRFS{kk,:}=irfs;
   FEVD{kk,:}=fevd;
   kk=kk+1
end



