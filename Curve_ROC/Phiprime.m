function result=Phiprime(X,t,Y,s)
if ((t>s)&&(X>Y))|| ((s>t)&&(Y>X)) 
    result=1;
elseif (t==s)|| (X==Y)
    result=0.5;
else
    result=0;
end