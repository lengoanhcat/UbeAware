function A = AUC(xs, ys)

% area under the curve defined by the _sequence_ of points [xs(i); ys(i)].
% xs can be either in increasing or in decreasing order; 
%

% force points to be in increasing order of x coordinate
if xs(end) < xs(1)
  xs = reverse(xs);
  ys = reverse(ys);
end

% compute area
A = 0;
for i = 2:length(xs)
  base = xs(i)-xs(i-1);
  height = (ys(i)+ys(i-1))/2;
  A = A + base*height;
end
