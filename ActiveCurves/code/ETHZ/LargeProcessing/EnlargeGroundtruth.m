function I = EnlargeGroundtruth(I, rf)

% Resizes groundtruth bbs I(ix).gtBB by scaling by rf (their center stays unchanged)
%
% Input:
% - I.gtBB(ix,:) = [minx miny maxx maxy]
% - rf = a scalar resize factor
%

for ix = 1:length(I)
  gtBB = I(ix).gtBB;
  for n = 1:size(gtBB,1)
    ctr = mean(reshape(gtBB(n,:),2,2)');
    w = gtBB(n,3)-gtBB(n,1);
    h = gtBB(n,4)-gtBB(n,2);
    I(ix).gtBB(n,:) = [ctr(1)-rf*w/2 ctr(2)-rf*h/2 ctr(1)+rf*w/2 ctr(2)+rf*h/2];
  end
end
disp(['Ground-truth BBs enlarged of factor ' num2str(rf)]);
