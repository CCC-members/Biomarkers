function AUCr = roc_stability_check(labels, scores, posclass, percent_in, ntimes, nx, xpos, grouplab, med, wrt_tit)
%labels: clase real de cada sujeto
%scores: clase estimada
%posclass: positive class
%percent_in: porciento de la muestra a tomar para calcular la ROC
%ntimes: # de veces que se va a realizar el procedim
%nx: #numero de subplots
%xpos: en cual subplot se va a plotear
%grouplab: etiqueta de los grupos. Es un cell
%med: es una etiqueta que se pone en el eje Y. El nombre de la medida
%wrt_tit: si es 1 se imprime un titulo en el subplot

plotdato = 0;

nclases = length(unique(labels));
if nclases == 2
    ny = 4;
else
    pair_ind = find(triu(ones(nclases),1));
    ny = length(pair_ind)+1;
end

if plotdato
    ny = ny+1;
    xpos = (xpos-1)*ny+1;
    subplot(nx, ny, xpos);
    plot(labels, scores, '.'); axis([min(labels)-0.1 max(labels)+0.1 min(scores)-0.2 max(scores)+0.2]);
    set(gca, 'xtick', unique(labels), 'xticklabel', grouplab); ylabel(med)
    title(['ROC Stability (' num2str(ntimes) ' Rep)'])
    pause(0.01)
    xpos = xpos+1;
end

if nclases == 2
    subplot(nx, ny, xpos); axis([0 1 0 1])
    [Xi,Yi, T, au] = perfcurve(labels, scores, posclass);
    plot(Xi,Yi); if wrt_tit, title('ROC whole sample'); end
    ind=find(Xi<=0.1); r = Yi(ind(end));
    text(0.3, 0.4, ['AUC 0.1 FP = ' num2str(r,'%3.2f')],'fontsize', 8)
    ind=find(Xi<=0.2); r = Yi(ind(end));
    text(0.3, 0.3, ['AUC 0.2 FP = ' num2str(r,'%3.2f')],'fontsize', 8)
    text(0.3, 0.2, ['AUC = ' num2str(au,'%3.2f')],'fontsize', 8)
    pause(0.01)
end

[X, Y, AUC] = resample_roc(labels, scores, posclass, percent_in, ntimes);

save mylastreasample
dbstop if error

if nclases == 2
    AUCr.AUC30 = zeros(2,1); AUCr.AUC50 = zeros(2,1);
    subplot(nx, ny, xpos+1); [AUCr.AUCt30, AUCr.AUCt50] = oneplot([AUC.total], wrt_tit, ['ROC Stab.(FP-100%)']);
    subplot(nx, ny, xpos+2); [AUCr.AUC30(1), AUCr.AUC50(1)] = oneplot([AUC.AUC10], wrt_tit, ['ROC Stab.(FP-10%)']);
    subplot(nx, ny, xpos+3); [AUCr.AUC30(2), AUCr.AUC50(2)] = oneplot([AUC.AUC20], wrt_tit, ['ROC Stab.(FP-20%)']);
    pause(0.01)
 
else
    subplot(nx, ny, xpos); [AUCr.AUCt30, AUCr.AUCt50] = oneplot([AUC.total], wrt_tit, ['ROC Stab.(AUC Tot)']);
    V = [AUC.pairwise];
    [ii, ij] = ind2sub([nclases, nclases], pair_ind);
    AUCr.AUC30 = zeros(size(V,1),1); AUCr.AUC50 = zeros(size(V,1),1);
    AUCr.pairs = cell(size(V,1), 1);
    for k=1:size(V,1)
        subplot(nx, ny, xpos+k); [AUCr.AUC30(k), AUCr.AUC50(k)]=oneplot(V(k,:), wrt_tit, ['ROC Stab.(' num2str(ii(k)) '-' num2str(ij(k)) ')']);
        AUCr.pairs{k} = ['(' num2str(ii(k)) '-' num2str(ij(k)) ')'];
    end
end


function [AUC30, AUC50]=oneplot(AUC, wrt_tit, tit)
vu = unique(AUC);
if length(vu) == 1
    Y1 = linspace(0.99,1, 100); X1 = ones(1,100).*vu.*100;
else
    [X1, Y1]=hist(AUC,length(vu));
    X1 = (length(AUC)-cumsum(X1))./length(AUC)*100;
end
bar(Y1,X1); hold on
ymin = min(Y1); ymax = max(Y1);
gap = (ymax - ymin)*0.2; if gap < eps; gap = 0.2; end
ymin = max(0, ymin - gap); ymax = ymax + gap;
axis([ymin ymax 0 100])
pc = 30; line(get(gca,'xlim'), [pc pc])
ii = find(X1>=pc);
xa = Y1(ii(end));
AUC30 = xa;
line([xa xa], [0 pc+10]); hh=text(xa, pc+10, num2str(xa, '%4.2f'), 'backgr', [1 1 0]);
pc = 50; line(get(gca,'xlim'), [pc pc])
ii = find(X1>=pc);
xa = Y1(ii(end));
AUC50 = xa;
line([xa xa], [0 pc+10]); text(xa, pc+10, num2str(xa, '%4.2f'), 'backgr', [1 1 0]);
if wrt_tit, title(tit); end
xa =get(gca, 'xlim');
xticks = linspace(xa(1), xa(2), 4);
set(gca, 'xtick', xticks, 'xticklabel', round(xticks.*1000)./1000, 'ytick', [0 30 50 70 100], 'ylim', [0 100]);
ylabel('Frequency (%)'); xlabel('AUC');
