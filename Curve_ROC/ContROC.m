function [OUTPUT] = ContROC(gldstd, test1, test2)
% This fuction is a traslated too R fuction (contROC.R) from Package ‘nonbinROC’
% Performs ROC-type analysis for continuous scale gold standard
% INPUT:
%      gldstd = vector giving the gold standard
%       test1 = vector giving a diagnostic test
%       test2 = optional vector giving another diagnostic test
% OUTPUT
%    OUTPUT.ROC = ROC curve
%    OUTPUT.Estimates = accuracy of the diagnostic test.
%    OUTPUT.StE = Standard Errors
%    OUTPUT.cov = covariance matriz between test1 and test2
%  
% References:
%   Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests’ 
%   accuracy when the gold standard is not binary. 
%   Academic Radiology, 12, 1198-1204.
% ----------------------------------------------------
% Author(s)
%   Paul Nguyen
%   April 17, 2009
% Traslated to Matlab by 
%   Msc Mayrim Vega Hernández
%   June, 2010.

%% Compute an estimator diagnostic test1 accuracy
N = length(gldstd); % total number of patiens in the study sample
vcomp1 = zeros(1,N); % Structural components of variance estimator
for i=1:N
    vcomp1(i) = sum(gldstd(i) == gldstd(setdiff(1:N,i)) | test1(i) == test1(setdiff(1:N,i)))/2 + ...
        sum(gldstd(i) < gldstd(setdiff(1:N,i)) & test1(i) < test1(setdiff(1:N,i))) + ...
        sum(gldstd(i) > gldstd(setdiff(1:N,i)) & test1(i) > test1(setdiff(1:N,i)));
end % sume case fi=1/2 if Xit = Xjs or t = s 
    % sume case fi=1   if Xit < Xjs or t < s 
    % sume case fi=1   if Xit > Xjs or t > s 
vcomp1 = vcomp1/(N - 1); % ROC; structural components of variance of fi^
acc1 = sum(vcomp1)/N; % Estimator of diagnostic test accuracy  
cov11 = sum((vcomp1 - acc1).^2)/(N/2 * (N/2 - 1)); %Estimate of variance for test1

%% compute an estimator of diagnostic test2 accuracy
if ~exist('test2')
    acc1 = max(acc1, 1 - acc1);
    
    OUTPUT.ROC = sort(vcomp1);
    OUTPUT.AreaEstimates = acc1;
    OUTPUT.StE = sqrt(cov11);
    OUTPUT.Cov = cov11;
else
    vcomp2 = zeros(1,N);
    for i=1:N
        vcomp2(i) = sum(gldstd(i) == gldstd(setdiff(1:N,i)) | test2(i) == test2(setdiff(1:N,i)))/2 + ...
            sum(gldstd(i) < gldstd(setdiff(1:N,i)) & test2(i) < test2(setdiff(1:N,i))) + ...
            sum(gldstd(i) > gldstd(setdiff(1:N,i)) & test2(i) > test2(setdiff(1:N,i)));
    end
        vcomp2 = vcomp2/(N - 1);
        acc2 = sum(vcomp2)/N;
        cov22 = sum((vcomp2 - acc2).^2)/(N/2 * (N/2 - 1));
        cov12 = sum((vcomp1 - acc1).* (vcomp2 - acc2))/(N/2 *(N/2 - 1));
        acc1 = max(acc1, 1 - acc1);
        acc2 = max(acc2, 1 - acc2);
        cov12 = abs(cov12);
             
        OUTPUT.ROC = [sort(vcomp1); sort(vcomp2)];
        OUTPUT.AreaEstimates = [acc1 acc2];
        OUTPUT.StE = [sqrt(cov11) sqrt(cov22)];
        OUTPUT.Cov = [cov11 cov12; cov12 cov22];
end
  