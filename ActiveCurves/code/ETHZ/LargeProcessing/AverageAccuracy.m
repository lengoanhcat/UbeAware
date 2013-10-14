function avg_acc = AverageAccuracy(scores, accs, threshs)

% returns acc(t) = average accuracy of detections with score >= thresh(i)
%

avg_acc = zeros(1,length(threshs));
for tix = 1:length(threshs)
  t = threshs(tix);
  scixs = (scores > t);
  if not(isempty(scixs))
    avg_acc(tix) = mean(accs(scixs));
  else
    avg_acc(tix) = 0;       % no average accuracy when no det at all ;)
  end
end
