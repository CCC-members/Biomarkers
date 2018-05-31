function [fit, type, alpha] = glmnet_regressors(datos, labels, alphas)
%ns: numero de sujetos
%datos: matriz de ns x nvar

ns = size(datos,1);
family = 'gaussian'; %'binomial'; %'gaussian';
type = 'response'; %'class'; %'response';

tt = 0;
%     figure
for alfa=alphas
    options=glmnetSet;
    options.alpha = alfa;
    tt=tt+1;
    %   subplot(3,3,tt)
    hh=figure;
    CVerr=cvglmnet(datos,labels,100,[],type, family,options,0);
    close(hh)
    mingcv(tt) = min(CVerr.cvm);
    minlambda(tt) = CVerr.lambda_min;
    disp(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
    %         title(['Alpha = ' num2str(alfa) '   GCVmin = '  num2str(min(CVerr.cvm))])
    %         drawnow
end

if length(alphas) > 1
    figure(2); delete(gca); plot(alphas, mingcv); pause(0.01)
end

% Con el alpha optimo, estimar el lambda optimo
[lambdamin,pos_alphamin]= min(mingcv);
alpha = alphas(pos_alphamin);
lmin = lambdamin;

% Con el lambda y el alpha escogidos, volver a calcular los betas
options=glmnetSet;
options.alpha = alpha;
options.nlambda = 1;
options.lambda = minlambda(pos_alphamin);

fit=glmnet(datos,labels, family, options);
