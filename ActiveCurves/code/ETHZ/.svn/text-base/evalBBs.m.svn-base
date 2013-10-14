function [res score]=evalBBs(BBs,testImg,pad)
% the padding should be bb pad, not image pad
trueBB=[];
trueScore=[];
falseBB=[];
falseScore=[];
% if it is negtive image
 score = zeros(length(BBs),1);
    for iBB = 1:length(BBs)
        score(iBB) = BBs(iBB).score;
    end
if ~isfield(testImg,'anno')
    res = zeros(length(BBs),1);
   
    return;
end
if isempty(testImg.anno)
    res = zeros(length(BBs),1);
   
    return;
end

anno  = testImg.anno;
anno(:,3) = testImg.anno(:,3)-testImg.anno(:,1);
anno(:,4) = testImg.anno(:,4)-testImg.anno(:,2);


predict = zeros(length(BBs),4);
for iBB = 1:length(BBs)
    predict(iBB,1) = (BBs(iBB).left+pad)/ BBs(iBB).scale;
    predict(iBB,2) = (BBs(iBB).top+pad)/ BBs(iBB).scale/BBs(iBB).asp;
    predict(iBB,3) = (BBs(iBB).wy-2*pad)/BBs(iBB).scale;
    predict(iBB,4) = (BBs(iBB).wx-2*pad)/BBs(iBB).scale/BBs(iBB).asp;
end
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
for iObj = 1:size(anno,1)
    ind = find(pascal(:,iObj)==1,1);
    pascal(ind+1:end,iObj) = 0;
end
res = sum(pascal,2)>0;

% for iBB = 1:length(BBs)
%     if( res(iBB) ==1)
%         BB = BBs(iBB);
%         BB.imgName = testImg.imgName;
%         trueBB = [trueBB;BBs(iBB)];
%         trueScore = [trueScore,BB.score];
%     else
%         BB = BBs(iBB);
%         BB.imgName = testImg.imgName;
%         falseBB = [falseBB;BBs(iBB)];
%         falseScore = [falseScore,BB.score];
%     end
% end
% BB = [trueBB;falseBB];
% score = [trueScore,falseScore];
end

