function result=theta(Xt,Xs)
result=0;
nt=length(Xt);
ns=length(Xs);
for t=1:nt
    for s=1:ns
        result=result+Phi(Xt(t),Xs(s));
    end
end
result=result/(nt*ns);