function [X,Y, AUC] = resample_roc(labels, score, posclass, percent_in, ntimes)

try
    delete('stopit.mat')
end

classes = unique(labels);
nclases = length(classes);
score = score(:);
labels = labels(:);
for k=1:nclases
    ind{k} = find(labels == classes(k));
    n(k) = length(ind{k});
    nsamp(k) = round(n(k)*percent_in/100);
end
if nclases==2
    X = 0:0.02:1;
    Y = zeros(length(X), ntimes);
    AUC = struct('total', cell(ntimes,1), 'AUC10', cell(ntimes,1), 'AUC20', cell(ntimes,1));
    i10 = find(X == 0.1);
    i20 = find(X == 0.2);
    for k=1:ntimes
        fprintf('.');
        if rem(k,80)==0,
            fprintf('%s', char(13));
        end
        i1 = randsample(n(1), nsamp(1));
        i2 = randsample(n(2), nsamp(2));
        [Xi,Yi, T, au] = perfcurve([labels(ind{1}(i1)); labels(ind{2}(i2))], [score(ind{1}(i1)); score(ind{2}(i2))], posclass, 'Xvals', X);
        Y(:, k) = aprox_roc(Xi,Yi,X);
        AUC(k).total = au;
        AUC(k).AUC10 = Y(i10,k);
        AUC(k).AUC20 = Y(i20,k);
        %     plot(Xi,Yi,'.', X, Y(:,k)); pause(0.1)

        %truco para detenerlo
        try
            load('stopit.mat')
            if stopit
                AUC = AUC(1:k);
                Y = Y(:, 1:k);
                break
            end
        end
        save resampling AUC X Y
    end
else
    X = []; Y = [];
    if nclases == 3
        L = [0 .5 .9 ; .5  0 .5 ; .9 .5 0];
    elseif nclases == 5
        L = [0 .5 .7 .8 .9;...
            .5  0 .5 .5 .5;...
            .7 .5  0 .5 .5;...
            .8 .5 .5  0 .5;...
            .9 .5 .5 .5  0];
    end
    tic
    for k=1:ntimes
        fprintf('.');
        if rem(k,80)==0,
            fprintf('%s', char(13));
        end
        xw = cell(nclases,1);
        ns = zeros(nclases,1);
        for h=1:nclases
            fit = score(ind{h});
            indh = randsample(n(h), nsamp(h));
            ns(h) = length(indh);
            xw{h}=fit(indh);
        end
        
        lsup = [0; cumsum(ns)];
        [AUC(k).total,pairwise]=ordinalROC(xw,L);
        AUC(k).pairwise = pairwise(find(triu(pairwise, 1)));
        
        %truco para detenerlo
        try
            load('stopit.mat')
            if stopit
                AUC = AUC(1:k);
                break
            end
        end
    end
    fprintf('%s', char(13));
    save resampling AUC X Y 
end
