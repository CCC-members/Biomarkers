function [selected, tested, elected, percent_selected, betas] = get_most_frequent(varsamples, fits, percent_selected, show_graph)
Nvar = max(varsamples(:));
Nsamples = size(varsamples,2);

elected = zeros(Nvar,1);
tested = zeros(Nvar,1);
betas = cell(Nvar,1);
for repet=1:Nsamples
    v = varsamples(:,repet);
    tested(v) = tested(v)+1;
    ii = find(fits{repet}.beta);
    v = v(ii);
    elected(v) = elected(v)+1;
    for h=1:length(v)
        betas{v(h)} = [betas{v(h)} fits{repet}.beta(ii(h))];
    end
end

if show_graph
    figure;
    subplot(221);hist(tested, Nvar);
    title(['Frequency of variables selected in ' num2str(Nsamples) ' samples']); xlabel('Number of Realizations'); ylabel('Number of Vars')
    subplot(222);hist(elected, Nvar)
    title('Number of times variables were selected'); xlabel('Number of Realizations Elected'); ylabel('Number of vars')
    elected = elected./tested*100;
    subplot(223); plot(elected, '.'); hold on
    line([1 Nvar], [50 50])
    line([1 Nvar], [40 40])
    title('Variables selection'); xlabel('Variables'); ylabel('Percent of times selected')
    subplot(224);hist(elected, Nvar)
    title('Percent of times variables were selected'); xlabel('Percent of Realizations selected'); ylabel('Number of vars')
end

selected = find(elected > percent_selected);
if isempty(selected),
    while isempty(selected) || (length(selected) < Nvar*0.03)
        percent_selected = percent_selected-1;
        selected = find(elected > percent_selected);
    end
end
