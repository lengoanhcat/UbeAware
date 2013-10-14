function [eer, th, auc] = PrecRecallPlot(corr_scores, false_scores, tot_corr, new_fig, draw_eer, col)

% Plots Precision/Recall curve and returns recall at the equal-error rate
%

% process arguments
if nargin < 3
  new_fig = true;
end
if nargin < 4
  draw_eer = false;
end
if nargin < 5
  col = [0 0 0];
end

[det_count false_count threshs] = ROCCount(corr_scores, false_scores);  

if new_fig
  figure;
end

ixs = det_count+false_count==0;
det_count(ixs) = 1;         % avoid div-by-0 in xs
xs = false_count./(det_count+false_count);
det_count(ixs) = 0;         % don't fake ys
ys = det_count/tot_corr;    % recall
xs = [1 xs];                % 1-precision
ys = [ys(1) ys];            % add leftmost point
threshs = [threshs(1) threshs];
plot(ys, 1-xs, 'LineWidth', 2, 'Color', col); hold on;
ylabel('Precision');
xlabel('Recall');
title('Precision/Recall plot');

if draw_eer
  plot([0 1], [1 0]);
end

% compute eer and auc
[trash eer_ix] =  min(abs(xs - (1-ys)));
eer = ys(eer_ix);
th = threshs(eer_ix);
auc = AUC(ys,1-xs);
