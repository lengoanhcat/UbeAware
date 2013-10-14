function [auc tp fp] = calc_auc(predicted,gt)
% <predicted> and <gt> should be two vectors, and their index correspond.
% <predicted> contain confidence values
% In <predicted>,value is larger for postive instances.
[so si] = sort(-predicted);
tp=cumsum(gt(si)>0)/sum(gt>0);
fp=cumsum(gt(si)<0)/sum(gt<0);
[uo,ui]=unique(so);
tp=[0;tp(ui);1];
fp=[0;fp(ui);1];

% compute lower envelope and area under curve
di=[true ; tp(2:end-1)~=tp(1:end-2) ; true];
x=fp(di);
y=tp(di);
auc=(x(2:end)-x(1:end-1))'*y(1:end-1);
%---------------------------------------
% CODE FROM PASCAL VOC2006 DEVELOP KIT
% (Modified)
%---------------------------------------