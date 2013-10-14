function [corr_dets, false_pos, corr_scores, false_scores, corr_overlaps] = CorrectDetections(C, gtBB, pascal)

% Partitions C into correct detections and false positives.
% A detection C(i) is correct if its bounding-box C(i).BB
% overlaps more than N% with one of the ground truth BBs gbBB(j,:) (bidirectional).
% At most one detection of C can be assigned to a ground truth BBs
% (multiple detections at the same location are counted as false-positives)
%
% if pascal == true -> uses the more severe PASCAL challenge criterion: intersection/union >= 0.5
%
% Output:
% corr_dets = indexes of correct detections of C
% false_pos = same for false positives
% corr_overlap = percentage area overlap with ground-truth BB (mean over two dirs)
%

if nargin < 3
  pascal = false;
end

% Convert all ground truth BBs to format of C.BB
% x_tl  x_br
% y_tl  y_br
%
n=0;
for bb = gtBB'
  n=n+1;
  gt(n).BB = [bb(1) bb(3); bb(2) bb(4)];
end

% Partitioning
corr_dets = [];
corr_scores = [];
corr_overlaps = [];
assigned_detections = zeros(1,length(C));  % detections already assigned to a ground truth BB
for j = 1:n          % loop over ground-truth BBs
  for i = find(not(assigned_detections))
    [ov1 a] = BBOverlap(C(i).BB,  gt(j).BB);
    [ov2 a] = BBOverlap(gt(j).BB, C(i).BB);
    corr = false;
    if pascal
      pcrit = a/(BBArea(C(i).BB)+BBArea(gt(j).BB)-a);  % intersection/union
      if pcrit >= 0.5  corr = true;  end
    else
      %if ov1 > 1/6 & ov2 > 1/6   % criterion until August 2005 report
      %if ov1 > 1/5 & ov2 > 1/5 % & ov1+ov2>2/3   % slightly stricter criterion used in ECCV06-subm
      if ov1 > 1/5 & ov2 > 1/5   % criterion for eccv06-final, bmvc06-subm, pami-kas sumbission
        corr = true;
      end
    end
    %
    if corr
      corr_dets = [corr_dets i];
      corr_scores = [corr_scores C(i).s];
      assigned_detections(i) = true;
      corr_overlaps = [corr_overlaps mean([ov1 ov2])];
      break;   % this gtBB is assigned
    end
  end
end

% false-positives data
false_pos = setdiff(1:length(C),corr_dets);
false_scores = [];
for i = false_pos
  false_scores = [false_scores C(i).s];
end
