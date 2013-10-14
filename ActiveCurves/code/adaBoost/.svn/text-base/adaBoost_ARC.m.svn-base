aucTable = zeros(nFold,nFeature);

matlabpool 4; % create workers for this task;

for iFold =1: nFold
pos_file_list = dir([cachePath '/' 'pos*.mat']);
num_pos = length(pos_file_list);
neg_file_list = dir([cachePath '/' 'neg*.mat']);
num_neg = length(neg_file_list);
num_data = num_pos + num_neg;

[train_p,test_p]=randSplit(num_pos,0.2);
[train_n,test_n]=randSplit(num_neg,0.2);
num_train_pos = length(train_p);
num_train_neg = length(train_n);
num_test_pos = length(test_p);
num_test_neg = length(test_n);
num_train = num_train_pos + num_train_neg;
num_test = num_test_pos + num_test_neg;
%% create data for active arcs
% get corresponding between feature dimensions and the mapping
load([cachePath '/' pos_file_list(1).name ]);
p1 = bound(1)+1; p2=bound(2);
p = p2-p1+1;

iPos = 1;
load([cachePath '/' pos_file_list(train_p(iPos)).name]);
crop_feature = fVec(p1:p2);
map  = find(crop_feature>-1e6);
p_map = p1:p2;
p_map = p_map(map);
p = length(p_map);


% allocate training data matrix
FullData= zeros(num_train,p,'single');

% alloate and fill the labels.
FullLabels = zeros(1,num_train);
FullLabels(1:num_train_pos)=1;
FullLabels(num_train_pos+1:end)=-1;

% fill out the training data matrix.
for iPos = 1:num_train_pos
    fprintf('train pos: %d of %d.\n ',iPos,num_train_pos);
    load([cachePath '/' pos_file_list(train_p(iPos)).name]);
    FullData(iPos,:) = fVec(p_map);
end
for iNeg = 1:num_train_neg
    fprintf('train neg %d of %d.\n',iNeg,num_train_neg)
    load([cachePath '/' neg_file_list(train_n(iNeg)).name]);
    FullData(num_train_pos+iNeg,:)=fVec(p_map);
end
save([outPath  '/ARCData_train_fold_' num2str(iFold) '.mat'],'train_p','train_n','FullData', 'FullLabels');

% train adaBoost model
iFeature = 0;
model =[] ;
data_weight = [];
[outClass,model,data_weight]=adaboost('train',FullData,FullLabels,model,data_weight,nFeature);

% clear data to save memory usage
FullData = [];
FullLabels = [];


FullTestLabels = zeros(1,num_test);
FullTestLabels(1:num_test_pos)=1;
FullTestLabels(num_test_pos+1:end)=-1;

outScoreMap = zeros(num_test,length(model));


save([outPath  '/ARCData_test_fold_' num2str(iFold) '.mat'],'test_p','test_n');
for iPos = 1:num_test_pos
	fprintf('test pos: %d of %d.\n', iPos, num_test);
	load([cachePath '/' pos_file_list(test_p(iPos)).name]);
	test_feature = fVec(p_map);
        for iFeature = 1:length(model)
            outScoreMap(iPos,iFeature) = adaboost('apply',test_feature,model(1:iFeature));
        end
end

for iNeg = 1:num_test_neg
       fprintf('test neg %d of %d.\n',iNeg + num_test_pos,num_test);
	load([cachePath '/' neg_file_list(test_n(iNeg)).name]);
	test_feature = fVec(p_map);
	for iFeature = 1:length(model)
	    outScoreMap(iNeg+num_test_pos,iFeature)=adaboost('apply',test_feature,model(1:iFeature));
	end
end

% calculate auc table
for iFeature=1: length(model);
  auc=calc_auc(outScoreMap(:,iFeature),FullTestLabels');
  aucTable(iFold,iFeature)=auc;
end

fileName = sprintf( '%s/ARC_Fold_%d_result.mat',outPath,iFold);
save(fileName,'aucTable','outScoreMap','model');
% visulize arc tempalte
vis_arc_template;
end

% plot the AUC curve;
meanAuc = mean(aucTable);
stdAuc = std(aucTable);
errorbar(1:nFeature,meanAuc,stdAuc,'r','linewidth',2)
axis([0 nFeature 0.5 1])

filename = sprintf('%s/ARC_adaBoost_nFeature_%d_nTrain_%d_nTest_%d_nFold_%d',outPath,nFeature,num_train_pos+num_train_neg,num_test_pos+num_test_neg,nFold);
title(filename)
xlabel('number  of features');
ylabel('AUC');
saveas(gcf,[filename '.fig']);
save([filename '.mat'],'aucTable','outScoreMap','model','p_map');

% close matlab workers
matlabpool close;
