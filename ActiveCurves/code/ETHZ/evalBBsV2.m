function [res score]=evalBBsV2(BBs,testImg)
% the padding should be bb pad, not image pad

%% rev 1, in this revision, bounding box is assumed to be object bounding box on input image coorindate system.
%% more specificly, BBs is a matrix, where each row is a bounding box of
%% [top left bottom right,... any intresting data ...., score of BB]
trueBB=[];
trueScore=[];
falseBB=[];
falseScore=[];

% if it is negtive image
nBB = size(BBs,1);
score = BBs(:,end);
if ~isfield(testImg,'anno')
    res = zeros(nBB,1);
    return;
end
if isempty(testImg.anno)
    res = zeros(nBB,1);
    return;
end



anno  = testImg.anno;
anno(:,3) = testImg.anno(:,3)-testImg.anno(:,1);
anno(:,4) = testImg.anno(:,4)-testImg.anno(:,2);

%[ left top width height]
predict = [BBs(:,2), BBs(:,1), BBs(:,4)-BBs(:,2), BBs(:,3)-BBs(:,1)];
intArea = rectint(predict,anno);

areaA = predict(:,3).*predict(:,4);
areaB = anno(:,3).*anno(:,4);
unionArea=zeros(size(intArea));
for i=1:length(areaB)
    unionArea(:,i)=areaA;
end
for i = 1:length(areaA)
    unionArea(i,:)=unionArea(i,:) + areaB';
end
unionArea = unionArea-intArea;

pascal = (intArea./unionArea)>0.5;
res = zeros(nBB,1)-1;
for iObj = 1:size(anno,1)
    ind = find(pascal(:,iObj)==1,1,'first');
    pascal(:,iObj) = 0;
    pascal(ind,:)=0;
    res(ind)=1;
end

res = res>0;
end

