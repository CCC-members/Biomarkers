function calling_glmnet_old_version(datos, grupo, edad, vnames, grouplab, parameters)

% Necessary parameter (with an example of possible default value) are:
% parameters.use_indfeat = 1;
% parameters.percent_var = 30;
% parameters.percent_subj = 70;
% parameters.percent_selected = 50;
% parameters.Nsamples = 500;
% parameters.nrepet = 500;
% parameters.percent_train_sample = 70;
% parameters.show_graphs = 1;
% parameters.family = 'gaussian'; %'binomial'; %'gaussian';
% parameters.type = 'response'; %'class'; %'response';
% parameters.include_age = 1;
% parameters.makeit_smooth = 1;
% The structure can have any number of other optional parameters

addpath(genpath('H:\_VIP\MatlabTools\common'))
addpath('H:\_VIP\MatlabTools\Classification\Fancy_varsel')
addpath('h:\_VIP\MatlabTools\glmnet')
addpath('H:\_VIP\MatlabTools\Curve_ROC')
addpath('H:\_VIP\MatlabTools\Variable Selection\Indfeat')

pnames = fieldnames(parameters);
for k=1:length(pnames)
    eval([pnames{k} ' = parameters.' pnames{k} ';']);
end



s=char(datetime); s = strrep(s, ':','-');

fsave_samples = [s 'samples_trials'];
fsave_rocs = [s 'ROCs_trials'];
fsave_selected = [s 'SelectedSpec'];
fname_clasif = [s 'clasificadores'];
fparameters = [s '_parameters.mat'];
fname_stable_roc_jpg = [s '-ROC_stable'];
fname_roc_jpg = [s '-ROC'];

if use_indfeat
    fsave_samples = [fsave_samples '_indfeat'];
    fsave_rocs = [fsave_rocs '_indfeat'];
    fsave_selected = [fsave_selected '_indfeat'];
    fname_clasif = [fname_clasif '_indfeat'];
end

if makeit_smooth
    fsave_samples = [fsave_samples '_fsmooth'];
    fsave_rocs = [fsave_rocs '_fsmooth'];
    fsave_selected = [fsave_selected '_fsmooth'];
    fname_clasif = [fname_clasif '_fsmooth'];
end
    
if include_age
    fsave_samples = [fsave_samples '_age'];
    fsave_rocs = [fsave_rocs '_age'];
    fsave_selected = [fsave_selected '_age'];
    fname_clasif = [fname_clasif '_age'];
end

fsave_samples = [fsave_samples '.mat'];
fsave_rocs = [fsave_rocs '.mat'];
fsave_selected = [fsave_selected '.mat'];
fname_clasif = [fname_clasif '.xls'];

outv = []; minlev = [];
if use_indfeat
    if makeit_smooth
        minlev = 1;
    else
        minlev = 1;
    end
    %Eliminar variables no importantes con indfeat
    labels = grupo;
    u = unique(labels);
    figure; p = 1; np = (length(u)*(length(u)-1))/2;
    outv = 1:size(datos,2);
    for k=1:length(u)-1
        ki = find(labels == u(k));
        for h=k+1:length(u)
            hi = find(labels == u(h));
            Xif = datos([ki(:); hi(:)],:);
            Yif = labels([ki(:); hi(:)]);
            Sig(p,:) = indfeat(Xif,Yif);
            subplot(np+1,1,p); plot(Sig(p,:));
            out_if = find(Sig(p,:) < minlev);
            outv = intersect(outv, out_if);
            p = p+1;
        end
    end
    subplot(np+1,1,np+1); plot(1:size(Sig,2), max(Sig,[],1), [1 size(Sig,2)], [minlev minlev]);
    
    %Quitar estas variables del analisis porque aportan poco
    datos(:,outv) = [];
end
pause(0.01);

save(fparameters, 'parameters', 'outv', 'minlev');

diagn = grupo';
alphas = linspace(0.01,0.9,10);  %1 es lasso
labels = diagn';
u = unique(labels);
NSubjects = zeros(length(u),1);
for k=1:length(u)
    NSubjects(k) = length(find(diagn == u(k)));
end
severidad =u;
Csize = NSubjects;
Nvar = size(datos,2);
datoscopy = datos;

[varsamples, subjsamples, fits, alpha, AUC] = find_biomarkers(...
    datos, edad, labels, NSubjects, alphas, Nvar, percent_train_sample, Nsamples, percent_var, percent_subj, show_graphs);
[selected, tested, elected, perc, betasfit] = get_most_frequent(varsamples, fits, percent_selected, show_graphs);

ii = find(round(elected(selected)) < 40);
selected(ii) = [];

[spsel, so] = sort(elected(selected), 'descend');
    
save(fsave_samples, 'varsamples', 'subjsamples','selected','vnames');
save(fsave_rocs, 'fits', 'type', 'alpha', 'AUC', 'Nvar', 'Nsamples')

save(fsave_selected, 'selected', 'elected');

datos = datos(:, selected);
%
fg = figure; maximize;
[fitsf, type, alphaf, AUCf, labelSubj, fitSubj] = ...
    aplicar_glmnet_resample(datos, labels, Csize, severidad, alphas, percent_train_sample, nrepet, edad,fg);
title('ROC classification by GLMNet resampling (100 repetitions) for the selected variables in the first step')
pause(0.01);

clear betas
figure; maximize
for k=1:length(fitsf)
    betas(:,k) = fitsf(k).beta;
    plot(fitsf(k).beta,'.'); hold on
end
plot(1:size(betas,1),betas', '.')
title(['Beta values of GLMNet for the ' num2str(size(betas,1)) ' vars of the first setp']); xlabel('Variables'); ylabel('Betas')
ind = find(betas == 0);
betas(ind) = NaN;
mbetas = nanmean(betas, 2);
pause(0.01);

if use_indfeat
    vars = setdiff(1:Nvar+length(outv), outv);
    vselected = vars(selected);

    st = {};
    for k=1:length(selected)
        st{k,1} = vnames{vselected(k)}{1};
        st{k,2} = num2str(frange(vnames{vselected(k)}{2}),'%5.2f');
        st{k,3} = num2str(spsel(k),'%5.2f');
        st{k,4} = num2str(mbetas(so(k)),'%10.6f');
    end
    pause(1);
    xlswrite(fname_clasif, st);

else
    vselected = selected;
    st = {};
    for k=1:length(selected)
        st{k,1} = vnames{vselected(k)};
        st{k,2} = num2str(spsel(k),'%5.2f');
        st{k,4} = num2str(mbetas(so(k)),'%10.6f');
    end
    pause(1);
    xlswrite(fname_clasif, st);
end


fit=fitsf(1);
fit.a0 = mean([fitsf.a0]);
fit.lambda = mean([fitsf.lambda]);
fit.beta = mbetas;
% fit.beta = mean([fitsf.beta],2);


save(fsave_selected, 'fitsf', 'selected', 'elected', 'fits', 'betas', 'fit', 'labelSubj', 'fitSubj');

Ns = size(datos,1);
myfit = zeros(Ns,1);
for k=1:Ns
    myfit(k)=glmnetPredict(fit,type,datos(k,:)); % make predictions
end


posclass = max(labels);
roc_plot(labels, myfit, posclass);
export_fig(fname_roc_jpg, '-jpg', '-nocrop')

%%
figure; maximize
roc_stability_check(labelSubj, fitSubj, posclass, 50, nrepet, 1, 1, grouplab, 'GLMNet Stable', 1)
export_fig(fname_stable_roc_jpg, '-jpg', '-nocrop')

