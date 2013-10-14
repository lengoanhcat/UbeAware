function sel_cv = RemoveDuplicates(cv)

% remove duplicate y values for the same x value
% on curve cv = [x1 x2 ... ; y1 y2 ...]
% by keeping the one with the highest y
%
% cv(1,:) must be sorted (in any direction)
%

cur_x = -inf;
sel_cv = [];
for i = 1:size(cv,2)
  if cv(1,i) ~= cur_x
    cur_x = cv(1,i);
    sel_cv = [sel_cv cv(:,i)];
  else
    sel_cv(:,end) = [cur_x; max(sel_cv(2,end),cv(2,i))];
  end
end
