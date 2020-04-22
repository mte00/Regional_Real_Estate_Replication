function [ir,wold]=impulse(A,CI,ndet,p,nsteps)

n=size(A,1);
A=varcompanion(A,ndet,n,p);    %A is without deterministic variables
J=[eye(n) zeros(n,(p-1)*n)];
ir=[];
wold=[];
for h=0:nsteps
    wold=cat(3,wold,(J*(A^h)*J')); 
    ir=cat(3,ir,(J*(A^h)*J')*CI);  %cat concatenates arrays: cat(dimension,A,B)
end

%multiplier=cumsum(wold,3);    %gives the interim multipliers