close all; clear
addpath('E:\willdesktop\COLOMBIA data\biomarkers\glmnet')
addpath('E:\willdesktop\COLOMBIA data\biomarkers\Curve_ROC')
addpath('E:\willdesktop\COLOMBIA data\biomarkers\Indfeat')

load region_2group
% load example
%%%EXPECTED VARIABLES%%
%- labels: a vector containing the subject's classification (1, 2, 3, etc), sorted from Normal to Pathologic or from Less Severe to the More Severe labels. Always starts in 1
%- data: a matrix of NSubjects x NVars.
%- age: is a vector with the ages of the subjects. It does not need to be the age. It can be any variable that distinguished the subjects.
%       It is used for the process of randomly selecting the subjects in the algorithms. It guarantees that in each repetition, in all groups the
%       subjects are balanced by age. If this condition is not needed, just create a vector of ONES: age = ones(Nsubjects, 1);

%%SETINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelsgroups = {'G1'; 'G2'};
% labelsgroups = {'G1'; 'G2'; 'G3'};
use_indfeat = 1; %perform pre-processing step to discard non-informative variables. It makes dimension reduction
alphas = linspace(0.01,0.9,10);  %elastic-net Regularization parameter: 0 is ridge; 1 is lasso
percent_var = 30;
% percent_var = 70;
percent_subj = 70;
percent_selected = 50;
Nsamples = 500;  %Change to 500 or 1000
nrepet = 500;    %Change to 500 or 1000
percent_train_sample = 70;
show_graphs = 1;
family = 'gaussian'; %'binomial'; %'gaussian';
type = 'response'; %'class'; %'response';

fsave_samples = 'random_samples.mat';
fsave_selected = 'selected_vars.mat';
fsave_rocs = 'robust_roc.mat';
fname_clasif = 'clasif.xlsx';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = unique(labels);
if use_indfeat %dimension reduction with indfeat
    figure; p = 1; np = (length(u)*(length(u)-1))/2;
    outv = 1:size(data,2);
    clear Sig
    for k=1:length(u)-1
        ki = find(labels == u(k));
        for h=k+1:length(u)
            hi = find(labels == u(h));
            Xif = data([ki(:); hi(:)],:);
            Yif = labels([ki(:); hi(:)]);
            Sig(p,:) = indfeat(Xif,Yif);
            subplot(np+1,1,p); plot(Sig(p,:));
%             out_if = find(Sig(p,:) < 1.2);
            out_if = find(Sig(p,:) < 1.5);
            outv = intersect(outv, out_if);
            p = p+1;
        end
    end
    subplot(np+1,1,np+1); plot(1:size(Sig,2), max(Sig), [1 size(Sig,2)], [1 1]);
    
    %Removing non-informative variables
    vars_in = setdiff(1:size(data,2), outv);
    disp(['INDFEAT: ' num2str(length(outv)) ' variables discarded out of ' num2str(size(data,2))]);
    data(:,outv) = [];

end


NSubjects = zeros(length(u),1);
for k=1:length(u)
    NSubjects(k) = length(find(labels == u(k)));
end
severity =u;
Csize = NSubjects;
Nvar = size(data,2);

[varsamples, subjsamples, fits, alpha, AUC] = find_biomarkers(...
    data, age, labels, NSubjects, alphas, Nvar, percent_train_sample, Nsamples, percent_var, percent_subj, show_graphs,use_indfeat);
[selected, tested, elected, perc, betasfit] = get_most_frequent(varsamples, fits, percent_selected, show_graphs);


ii = find(round(elected(selected)) < percent_selected);
selected(ii) = [];

[spsel, so] = sort(elected(selected), 'descend');
    
save(fsave_samples, 'varsamples', 'subjsamples','vars_in', 'selected', 'elected');
save(fsave_rocs, 'fits', 'type', 'alpha', 'AUC', 'vars_in', 'Nvar', 'Nsamples')

save(fsave_selected, 'vars_in', 'selected', 'elected');

data = data(:, selected);
%
fg = figure;
[fitsf, type, alphaf, AUCf, labelSubj, fitSubj] = ...
    aplicar_glmnet_resample(data, labels, Csize, severity, alphas, percent_train_sample, nrepet, age,fg);
title('ROC classification by GLMNet resampling (100 repetitions) for the selected variables in the first step')

clear betas
figure
for k=1:length(fitsf)
    betas(:,k) = fitsf(k).beta;
    plot(fitsf(k).beta,'.'); hold on
end
plot(betas, '.')
title(['Beta values of GLMNet for the ' num2str(length(betas)) ' vars of the first setp']); xlabel('Variables'); ylabel('Betas')


if use_indfeat
    vselected = vars_in(selected);
else
    vselected = 1:Nvar;
end
xlswrite(fname_clasif, vselected);


%Getting the final classifier equation
fit=fitsf(1);
fit.a0 = mean([fitsf.a0]);
fit.lambda = mean([fitsf.lambda]);
fit.a0 = mean([fitsf.a0]);
fit.beta = mean([fitsf.beta],2);

save(fsave_selected, 'vars_in', 'selected', 'elected', 'fits', 'betas', 'fit', 'labelSubj', 'fitSubj');

Ns = size(data,1);
myfit = zeros(Ns,1);
for k=1:Ns
    myfit(k)=glmnetPredict(fit,type,data(k,:)); % make predictions
end


posclass = max(labels);
roc_plot(labels, myfit, posclass);

%%
figure
roc_stability_check(labelSubj, fitSubj, posclass, 40, 500, 1, 1, labelsgroups, 'GLMNet Stable', 1)
