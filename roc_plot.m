function [AUCc, AUCc10, AUCc20] = roc_plot(labels, myfit, posclass)
[X1,Y1, T, AUCc] = perfcurve(labels, myfit, posclass);
figure; set(gcf, 'color', 'w');
plot(X1,Y1);
xlabel('False positive rate', 'fontsize', 14); ylabel('True positive rate', 'fontsize', 14)
title('ROC for classification by GLMNet', 'fontsize', 14)
text(0.4, 0.5, ['AUC Total= ' num2str(AUCc,'%6.2f')], 'fontsize', 14)
ii = find(X1 <= 0.1);
ii = ii(end);
AUCc10 = Y1(ii);
text(0.4, 0.4, ['AUC 0.1% FP = ' num2str(AUCc10,'%6.2f')], 'fontsize', 14)
ii = find(X1 <= 0.2);
ii = ii(end);
AUCc20 = Y1(ii);
text(0.4, 0.3, ['AUC 0.2% FP = ' num2str(AUCc20,'%6.2f')], 'fontsize', 14)
