function [GTO, k] = SplitAndThinGTO(GTO, gto_type)

% splits multiple instances in GTO,
% and thins them to 1 pixel wide.
%
% Returns also number of instances k.
%
% Supports both gto_type 'vitto' and 'marcin'
%

% deal with multiple instances
if strcmp(gto_type,'vitto')
  GTO = SplitImage(GTO);
elseif strcmp(gto_type,'marcin')
  % nothing to do
else
  error([mfilename ': unknown gto type: ' gto_type]);
end
k = size(GTO,3);                                                       % number of instances

% thin GTO
for kix = 1:k
  GTO(:,:,kix) = bwmorph(GTO(:,:,kix),'thin',inf);
end
