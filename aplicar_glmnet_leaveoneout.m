function aplicar_glmnet_leaveoneout(datos, labels, Csize, severidad, alphas)
%ns: numero de sujetos
%datos: matriz de ns x nvar

%Se hace un ciclo para cada sujeto para estimar el alpha y el lambda
%optimos y evaluar al sujeto

ns = size(datos,1);
family = 'gaussian'; %'binomial'; %'gaussian';
type = 'response'; %'class'; %'response';

for k=1:ns
    disp(['Estimating subject ' num2str(k)])
    datok = datos(setdiff([1:ns],k),:);
    labelk = labels(setdiff([1:ns],k));

    tt = 0;
    %     figure
    for alfa=alphas
        options=glmnetSet;
        options.alpha = alfa;
        tt=tt+1;
        %   subplot(3,3,tt)
        hh=figure;
        CVerr=cvglmnet(datok,labelk,100,[],type, family,options,0);
        close(hh)
        mingcv(tt) = min(CVerr.cvm);
        minlambda(tt) = CVerr.lambda_min;
        disp(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
        %         title(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
        %         drawnow
    end
    
    if length(alphas) > 1
        figure(2); delete(gca); plot(alphas, mingcv); title(['Sujeto ' num2str(k)]); drawnow
    end
    
    % Con el alpha optimo, estimar el lambda optimo
    [lambdamin,pos_alphamin]= min(mingcv);
    alpha(k) = alphas(pos_alphamin);
    lmin(k) = lambdamin;
    
    % Con el lambda y el alpha escogidos, predecir el sujeto k
    options=glmnetSet;
    options.alpha = alpha(k);
    options.nlambda = 1;
    options.lambda = minlambda(pos_alphamin);
    
    disp(['Estimating subject ' num2str(k)])
    fit1=glmnet(datok,labelk, family, options);
    fits(k) = fit1;
    
    fitSubj(k)=glmnetPredict(fit1,type,datos(k,:)); % make predictions
    
end

%% ROC
%esto es para crear la variable que le entra a ordinalROC, que tiene que ser un arreglo de celulas,
%donde cada celula es un grupo: Ausencia de Lesion, lesion Clasica, Lesion No Clasica

if length(Csize) == 2
    L = [0 .9; .9 0];
else
    L = [0 .5 .9; .5  0 .5 ; .9 .5 0];
end

lsup = [0; cumsum(Csize)];

xw = cell(length(Csize),1);

for k=1:length(Csize)
    ind=find(labels==severidad(k));
    xw{k}=fitSubj(ind);
end

[thetapp,pairwise]=ordinalROC(xw,L)


%%

[X1,Y1] = perfcurve(labels, fitSubj, 2);
plot(X1,Y1);
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC for classification by GLMNet')
text(0.2, 0.5, ['AUC Total= ' num2str(pairwise(1,2),'%6.2f')])
ii = find(X1 <= 0.1);
ii = ii(end);
text(0.2, 0.4, ['AUC 0.1% FP = ' num2str(Y1(ii),'%6.2f')])
ii = find(X1 <= 0.2);
ii = ii(end);
text(0.2, 0.3, ['AUC 0.2% FP = ' num2str(Y1(ii),'%6.2f')])

% xx=[fits.beta];
% zz=sum(xx,2);
% yy =zeros(nd,nf);
% yy(quedan) = zz;
% % yy=reshape(yy,19,nf-length(quitar));
% figure; imagesc([1:nf]*freqres, 1:nd, yy); colorbar


