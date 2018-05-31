function [fits, type, alpha, AUC, labelSubj, fitSubj] = ...
    aplicar_glmnet_resample(datos, labels, Csize, severidad, alphas, percent_train_sample, nrepet, var_indep,fig)
%ns: numero de sujetos
%datos: matriz de ns x nvar

%Se selecciona una muestra de entrenamiento y una muestra de test al azar. Se estima el alpha y el
%lambda optimos usando la muestra de entrenamiento. Luego se evaluan los sujetos de la muestra
%de test. Se repite el proced nrepet veces

AUC.total = 0;
AUC.AUC10 = 0;
AUC.AUC20 = 0;

ns = size(datos,1);
family = 'gaussian'; %'binomial'; %'gaussian';
type = 'response'; %'class'; %'response';

severidad = unique(labels);
if length(severidad) > 4 %%== length(labels)
    severidad = []; %%ones(length(labels),1);
    lval = labels;
    ngroups = 1;
    NSubjects = length(labels);
    subjects{1} = 1:NSubjects;
else
    lval = unique(labels);
    ngroups = length(lval);
    for k=1:ngroups
        subjects{k} = find(labels == lval(k));
        NSubjects(k) = length(subjects{k});
    end
end

ns = zeros(ngroups,1);
ntrain = zeros(ngroups,1);
ntest = zeros(ngroups,1);
for k=1:ngroups
    ns(k) = NSubjects(k);
    ntrain(k) = round(ns(k).*percent_train_sample/100);
    ntest(k) = ns(k) - ntrain(k);
end

posfit = 1;
for k=1:nrepet
    disp(['Repeticion ' num2str(k)])
    
    %escoger muestra de train, tratar de que sea uniforme en el eje de las edades
    datotrain = []; labeltrain = [];
    datotest = []; labeltest = [];
    for h=1:ngroups
        ii = subjects{h};
        ind = randsample(ns(h), ntrain(h));
        [se, ie] = sort(var_indep(ii));
        itrain = ii(ie(ind));
        itest = ii(ie(setdiff(1:ns(h), ind)));
        
        datotrain = [datotrain; datos(itrain,:)];
        labeltrain = [labeltrain; labels(itrain)];
        datotest = [datotest; datos(itest,:)];
        labeltest = [labeltest; labels(itest)];
    end

    tt = 0;
    %     figure
    for alfa=alphas
        options=glmnetSet;
        options.alpha = alfa;
        tt=tt+1;
        %   subplot(3,3,tt)
        hh=figure;
        CVerr=cvglmnet(datotrain,labeltrain,100,[],type, family,options,0);
        close(hh)
        mingcv(tt) = min(CVerr.cvm);
        minlambda(tt) = CVerr.lambda_min;
        disp(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
        %         title(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
        %         drawnow
    end
    
%     if length(alphas) > 1
%         figure(20);
%         delete(gca); plot(alphas, mingcv); title(['Sujeto ' num2str(k)]); pause(0.01)
%     end
%     
    % Con el alpha optimo, estimar el lambda optimo
    [lambdamin,pos_alphamin]= min(mingcv);
    alpha(k) = alphas(pos_alphamin);
    lmin(k) = lambdamin;
    
    % Con el lambda y el alpha escogidos, volver a calcular los betas
    options=glmnetSet;
    options.alpha = alpha(k);
    options.nlambda = 1;
    options.lambda = minlambda(pos_alphamin);
    
    disp(['Estimating subject ' num2str(k)])
    fit1=glmnet(datotrain,labeltrain, family, options);
    fits(k) = fit1;
    if length(fit1) < 15
        disp(['Betas: ' num2str(fit1.beta')]);
    end
    
    %Predecir los sujetos que quedaron fuera
    fitSubj(posfit:posfit+size(datotest,1)-1)=glmnetPredict(fit1,type,datotest); % make predictions
    labelSubj(posfit:posfit+size(datotest,1)-1) = labeltest;
    posfit = posfit + size(datotest,1);
    
end

%% ROC

if ~isempty(severidad)
    %%%??????
    % indout = find(fitSubj == 0);
    % labelSubj(indout) = [];
    % fitSubj(indout) = [];
    
    for k = 1:length(severidad)
        Csize(k) = length(find(labelSubj == severidad(k)));
    end
    
    
    %esto es para crear la variable que le entra a ordinalROC, que tiene que ser un arreglo de celulas,
    %donde cada celula es un grupo: Ausencia de Lesion, lesion Clasica, Lesion No Clasica
    
    % L = [0 .9; .9 0];
    if length(Csize) > 2
        if length(Csize) == 3
            L = [0 .5 .9; .5  0 .5 ; .9 .5 0];
        else
            L = 0.9*ones(length(Csize));
            L(1:size(L,1)+1:end) = 0;
        end
        
        lsup = [0; cumsum(Csize)];
        
        xw = cell(length(Csize),1);
        
        for k=1:length(Csize)
            ind=find(labelSubj==severidad(k));
            xw{k}=fitSubj(ind);
        end
        
        [AUC.total,pairwise]=ordinalROC(xw,L);
        AUC.pairwise = pairwise(find(triu(pairwise, 1)));
        
    else
        
        if length(severidad) > 1
%             figure(19); clf;
            posclass = max(severidad);
            [X1,Y1, T, AUC.total] = perfcurve(labelSubj(:)+1-min(severidad), fitSubj(:)+1-min(severidad), posclass+1-min(severidad));
%             plot(X1,Y1); set(gca, 'fontsize', 12);
%             xlabel('False positive rate', 'fontsize', 12); ylabel('True positive rate', 'fontsize', 12)
%             title('ROC for classification by GLMNet', 'fontsize', 14)
%             text(0.2, 0.5, ['AUC Total= ' num2str(AUC.total,'%6.2f')], 'fontsize', 12)
            ii = find(X1 <= 0.1);
            ii = ii(end);
            AUC.AUC10 = Y1(ii);
%             text(0.2, 0.4, ['AUC 0.1% FP = ' num2str(AUC.AUC10,'%6.2f')], 'fontsize', 12)
            ii = find(X1 <= 0.2);
            ii = ii(end);
            AUC.AUC20 = Y1(ii);
%             text(0.2, 0.3, ['AUC 0.2% FP = ' num2str(AUC.AUC20,'%6.2f')], 'fontsize', 12)
        end
    end
    
    % xx=[fits.beta];
    % zz=sum(xx,2);
    % yy =zeros(nd,nf);
    % yy(quedan) = zz;
    % % yy=reshape(yy,19,nf-length(quitar));
    % figure; imagesc([1:nf]*freqres, 1:nd, yy); colorbar
    
else
end

