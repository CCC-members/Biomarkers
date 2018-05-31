function [varsamples, subjsamples, fits, alpha, AUC] = find_biomarkers(...
    datos, edad, labels, NSubjects, alphas, Nvar, percent_train_sample, Nsamples, percent_var, percent_subj, show_graphs, use_indfeat)

%Aqui la edad no se mete en los analisis, se usa a la hora de generar las muestras aleatorias, para que balanceadas por edad. Si no interesa, se
%puede poner de 1 a N u otra cosa

Nvar = size(datos,2);

severidad = unique(labels);
if length(severidad) > 4 %%== length(labels)
    severidad = []; %%ones(length(labels),1);
    NSubjects = sum(NSubjects);
end

varsamples = generate_samples(Nvar, percent_var, Nsamples);
Ngroups = length(NSubjects);
for k=1:Ngroups
    subjsamples{k} = generate_samples(NSubjects(k), percent_subj, Nsamples);
end

if show_graphs
    figure
    subplot(Ngroups+1,1,1); hist(varsamples(:), Nvar); title(['Variables selection ' num2str(Nvar) ' Variables.  ' num2str(Nsamples) ' samples'])
    for k=1:Ngroups
        subplot(Ngroups+1,1,k+1); hist(subjsamples{1}(:), NSubjects(1));
        title(['Subjects Sample ' num2str(k) ' selection ' num2str(NSubjects(k)) ' Subjects.  ' num2str(Nsamples) ' samples'])
    end
end

fits = cell(Nsamples, 1); alpha = zeros(Nsamples, 1);

nrepet = 1;

% fg = figure;
for repet=1:Nsamples
    fprintf('Sampling  %d of %d\n', repet, Nsamples);
    data = [];  Ns = zeros(length(NSubjects),1); labsamp = [];
    if isempty(severidad)
        ll = 1;
        labsamp = [labsamp; labels(subjsamples{ll}(:,repet))];
        Ns = size(subjsamples{ll},1);
        data = [data; datos(subjsamples{ll}(:,repet),varsamples(:,repet))];
    else
        for ll=1:length(severidad)
            ii = find(labels == severidad(ll));
            labsamp = [labsamp; labels(ii(subjsamples{ll}(:,repet)))];
            Ns(ll) = size(subjsamples{ll},1);
            data = [data; datos(ii(subjsamples{ll}(:,repet)),varsamples(:,repet))];
        end
    end
    
    if use_indfeat
        %Eliminar variables no importantes con indfeat
        u = unique(labsamp);
        p = 1;
        outv = 1:size(data,2);
        clear Sig
        for k=1:length(u)-1
            ki = find(labsamp == u(k));
            for h=k+1:length(u)
                hi = find(labsamp == u(h));
                Xif = data([ki(:); hi(:)],:);
                Yif = labsamp([ki(:); hi(:)]);
                Sig(p,:) = indfeat(Xif,Yif);
                out_if = find(Sig(p,:) < 1.2);
                outv = intersect(outv, out_if);
                p = p+1;
            end
        end
        %Quitar estas variables del analisis porque aportan poco
        disp(['-->> REPEATED INDFEAT: ' num2str(length(outv)) ' variables removed out of de ' num2str(size(data,2))])
        NV = size(data,2);
        quedaron = setdiff(1:NV, outv);
        data(:,outv) = [];
    end
    
    fg = [];
%     figure(fg)
    [fits{repet}, type, alpha(repet), AUC(repet)] = ...
        aplicar_glmnet_resample(data, labsamp, Ns, severidad, alphas, percent_train_sample, nrepet, edad, fg);
%     title(['ROC classification by GLMNET. Rep=' num2str(repet)]);
    betas = zeros(NV, 1);
    betas(quedaron) = fits{repet}.beta;
    fits{repet}.beta = betas;
%      
end

if (length(unique(labels)) == 2) & show_graphs
    figure;
    subplot(Ngroups+1,1,1); hist([AUC.total], 100); title(['Hist of AUC Total with all variables. ' num2str(Nsamples) ' resampling']);
    subplot(312); hist([AUC.AUC10], 100); title(['Hist of AUC 10% FP with all variables. ' num2str(Nsamples) ' resampling']);
    subplot(313); hist([AUC.AUC20],100); title(['Hist of AUC 20% FP with all variables. ' num2str(Nsamples) ' resampling']);
end
