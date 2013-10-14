function Sacc = AccumulateStats(S)

% combines two or more stat records S(k) into a single one,
% accumulating all stats.
%

if length(S) < 2
  Sacc = S;
  return;
end

Sacc = S(1);
for Scur = S(2:end)
  Sacc.corr_scores = [Sacc.corr_scores Scur.corr_scores];
  Sacc.false_scores = [Sacc.false_scores Scur.false_scores];
  if isfield(Scur,'corr_overlaps');
    Sacc.corr_overlaps = [Sacc.corr_overlaps Scur.corr_overlaps];
  end
  Sacc.pos_scores = [Sacc.pos_scores Scur.pos_scores];
  Sacc.neg_scores = [Sacc.neg_scores Scur.neg_scores];
  Sacc.tot_corr = Sacc.tot_corr + Scur.tot_corr;
  Sacc.tot_test_imgs = Sacc.tot_test_imgs + Scur.tot_test_imgs;
  Sacc.tot_pos = Sacc.tot_pos + Scur.tot_pos;
  Sacc.tot_neg = Sacc.tot_neg + Scur.tot_neg;
end
