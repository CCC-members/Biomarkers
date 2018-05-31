function [fits, fitSubj, labelSubj, mingcv, lambdamin, alpha, lmin] = glmnetgcv_onesample(datotrain, labeltrain, datotest, type, family, labeltest)
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
alpha = alphas(pos_alphamin);
lmin = lambdamin;

% Con el lambda y el alpha escogidos, volver a calcular los betas
options=glmnetSet;
options.alpha = alpha;
options.nlambda = 1;
options.lambda = minlambda(pos_alphamin);

fit1=glmnet(datotrain,labeltrain, family, options);
fits = fit1;
if length(fit1) < 15
    disp(['Betas: ' num2str(fit1.beta')]);
end

%Predecir los sujetos que quedaron fuera
fitSubj=glmnetPredict(fit1,type,datotest); % make predictions
labelSubj = labeltest;

