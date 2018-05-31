function Ya = aprox_roc(X1, Y1, Xa)
Ya = zeros(length(Xa),1);
for ir =1:length(Xa)
    ii = find(X1 <= Xa(ir));
%     [aa,ii] = min(abs(X1 - Xa(ir)));
    Ya(ir) = Y1(ii(end));
end
