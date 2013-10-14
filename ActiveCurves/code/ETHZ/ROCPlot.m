function [eer, auc] = ROCPlot(S, C, new_fig, draw_eer, col)

% Plots ROC curve of scores S,
% given that a correct classif is
% C(i)==1 --> S(i)>thresh.
% The threshold is varied over all S values.
%
% Returns equal error rate eer and area under curve auc
%

if nargin < 3
  new_fig = true;
end

if nargin < 4
  eer = false;
end

if nargin < 5
  col = [0 0 0];
end

%Smin = min(S);
%Smax = max(S);
%ts = Smin:(Smax-Smin)/1000:Smax;  % produces a visually smoother, but less accurate plot
ts = sort(S);
[det_rate, false_pos] = DetectionRates(S, C, ts);
if new_fig
  figure;
end

xs = false_pos;
ys = det_rate;
plot(xs, ys, 'LineWidth', 2, 'Color', col); hold on;

if draw_eer
  plot([0 1], [1 0]);
end
xlabel('Percentage false-positives');
ylabel('Percentage correct detections');
title('ROC plot');

[trash eer_ix] =  min(abs(xs - (1-ys)));
eer = ys(eer_ix);
auc = AUC(xs,ys);
