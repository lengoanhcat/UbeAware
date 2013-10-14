function [det_rates, false_poss] = DetectionRates(S, C, ts)

% detection rates and false positives
% for thresholds ts(i)

det_rates = [];
false_poss = [];

for t = ts
  [dr fp] = DetectionRate(S, C, t);
  det_rates = [det_rates dr];
  false_poss = [false_poss fp];
end
