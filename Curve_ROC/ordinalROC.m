function [thetapp,pairwise]=ordinalROC(X,L)
%X es cell array de T estandard oro ordenados
%L matriz costo TxT

T=length(X);
n=zeros(T,1);
pairwise=zeros(T);
for t=1:T;n(t)=length(X{t});end
w=n*n'./sum(sum(triu(n*n',1)));
thetapp=0;
for t=1:T
    for s=t+1:T
        pairwise(t,s)=theta(X{t},X{s});
        thetapp=thetapp+w(t,s)*L(t,s)*(1-pairwise(t,s));
    end
end
thetapp=1-thetapp;

