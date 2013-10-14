function [M S] = FuseCurves(curves, max_x, plot_curves)

% Fuse curves{i} = [x1 x2 x3 ...; y1 y2 y3 ...]
% into a mean curve M, with standard-devations S.
% Also plots the mean curve and the resampled input curves.
%
% This is an adaptation of Frederic Jurie's code.
%

% parse arguments
if nargin < 3
  plot_curves = false;
end

NC = length(curves);           % number of curves
N=30+1;                        % number of points sampled on each curve
val=linspace(0,max_x,N);
MR=zeros(NC,N);                % sampled points on each curve (now all sampled at the same x coords !)

for ic = 1:NC                  % loop over curves
  cv = curves{ic};             % org curve points
  cv = sort(cv',1)';           % must be sorted in ascending x order for this routine to work
  cv = RemoveDuplicates(cv);   % keep only highest y for each x
  cv = cv';
  np = size(cv);
  %
  for i = 1:N                  % loop over points to sample
    x = val(i);
    l1 = find(cv(:,1)<=x);
    l2 = find(cv(:,1)>=x);
    if size(l1,1) >0
      i1=max(l1);
      x1=cv(i1,1); v1=cv(i1,2);
    else
      x1=cv(1,1); v1=cv(1,2);
    end
    if size(l2,1) > 0
      i2=min(l2);
      x2=cv(i2,1); v2=cv(i2,2);
    else
      x2=cv(np,1); v2=cv(np,2);
    end;
    %
    d = x2-x1;
    if d~=0
     % v=(x-x1)/d*v2+(x2-x)/d*v1;   % mathematically correct line (thanks to Tingting); leads to tiny deviations from plots as in the papers)
     v=(x-x1)/d*v1+(x2-x)/d*v2;     % original line (leads to idential plots as in papers)
    else
      v=v1;
    end;
    MR(ic,i)=v;
  end % loop over sampled points
end % loop over curves

M=mean(MR);
S=std(MR);

if plot_curves
  figure; hold on;
  title('Resampled average curve with standard deviation');
  errorbar(val,M,S);
  axis([0 max_x 0 1]);
  grid on;

  figure;
  title('Resampled input curves');
  for ic=1:NC
    hold on
    plot(val,MR(ic,:), 'color', rand(1,3),'linewidth',2);
  end;
end

% format outputs;
M = [val; M];

