function [det_count, falsepos_count, thresholds] = ROCCount(pos, neg)

% comment to be written

thresholds = [sort([neg pos]) max([neg pos])];    % make a point at every change
det_count = length(pos);
falsepos_count = length(neg);
for t = thresholds
  det_count = [det_count sum(pos>t)];              % let it to '>', not '>=' --> higher plot
  falsepos_count = [falsepos_count sum(neg>t)];
end
thresholds = [-inf thresholds];
