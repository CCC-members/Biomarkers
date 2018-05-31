load Datos
%load datanew
%load mahalabien
%load datacomplete
%load data %sin mahalnobis
%load heart
nvar=size(data,2);
data(:,1)=data(:,1)+(1-min(data(:,1)));
niveles=length(unique(data(:,1)));
%L=[0 0.5 1;0.5 0 0.5;1 0.5 0];
desplaz=2;
L=toeplitz([0 0.5 1]);
%desplaz=0;
%desplaz=1;
thetapp=zeros(1,nvar-1);
etiqueta=cell(1,nvar);
for k=1:nvar;etiqueta{k}=textdata{1,k+desplaz};end
for k=2:nvar;
    X=cell(1,niveles);
    for t=1:niveles
        temp=data(t==(data(:,1)),k);
        temp(isnan(temp))=[];
        X{t}=temp;
     end
    [thetapp(k-1),pairwise]=ordinalROC(X,L);
    disp([ etiqueta{k},'=',num2str(thetapp(k-1))]);
   disp(pairwise)

end
[x,ind]=sort(thetapp,'descend');
barh(thetapp(ind));
temp=cell(1,nvar-1);
for k=1:nvar-1;temp{k}=etiqueta{ind(k)+1};end
set(gca,'YTickLabel',temp)
