function cids = ThresholdDetections(C, thresh)

% cid in cids iff C(cid).s > thresh
%

cids = [];
for cid = 1:length(C)
  if C(cid).s > thresh
    cids = [cids cid];
  end
end
