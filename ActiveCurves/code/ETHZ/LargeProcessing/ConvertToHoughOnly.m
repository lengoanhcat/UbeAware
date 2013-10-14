function NH = ConvertToHoughOnly(H, model)

% Converts H to plain hypotheses after the Hough stage,
% i.e. get the score from H(hix).basis(4),
% and the bb by transl+scale model by H(hix).basis(1:3)
%

if isempty(H)
  NH = [];
  return;
end

for hix = 1:length(H)
  basis = H(hix).basis;
  NH(hix).s = basis(4);
  vs(hix) = basis(4);
  ctr = mean(model.bb')';
  t = basis(1:2)'-ctr;
  s = basis(3);
  M_bb = model.bb;
  I_bb = [(M_bb(1,:)-ctr(1))*s + ctr(1) + t(1);          % model BB trafoed by t,sl
          (M_bb(2,:)-ctr(2))*s + ctr(2) + t(2)];
  NH(hix).bb = I_bb;
end

[trash order] = sort(-vs);
NH = NH(order);
