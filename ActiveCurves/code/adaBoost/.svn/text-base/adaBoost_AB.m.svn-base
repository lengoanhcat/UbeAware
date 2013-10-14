aucTable = zeros(nFold,nFeature);

matlabpool 4 % open paralle computation works, valid for matlab 2010a or later.
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

%% Create data matrix for active basis, the matrix is of size (np+nn) * p, where
%% np and nn are number of positive and negative examples, and p is dimension of 
%% feature.

% get bound for active basis feature;
	load([cachePath '/' pos_file_list(1).name ]);
	p1 = 0; p2=bound(1);
	p = p2-p1;

% Allocate the data matrix for training
FullData= zeros(num_train,p,'single');
FullLabels = zeros(1,num_train);
FullLabels(1:num_train_pos)=1;
FullLabels(num_train_pos+1:end)=-1;
% Fill out the training data matrix;
for iPos = 1:num_train_pos
    fprintf('train pos: %d of %d.\n ',iPos,num_train_pos);
    load([cachePath '/' pos_file_list(train_p(iPos)).name]);
    FullData(iPos,:) = fVec(p1+1:p2);
end
for iNeg = 1:num_train_neg
    fprintf('train neg %d of %d.\n',iNeg,num_train_neg)
    load([cachePath '/' neg_file_list(train_n(iNeg)).name]);
    FullData(num_train_pos+iNeg,:)=fVec(p1+1:p2);
end

% train the adaBoost model based using active basis as features
iFeature = 0;
model =[] ;
data_weight = [];
[outClass,model,data_weight]=adaboost('train',FullData,FullLabels,model,data_weight,nFeature);



% Allocate the data matrix for testing
FullTestData= zeros(num_test,p,'single');
FullTestLabels = zeros(1,num_test);
FullTestLabels(1:num_test_pos)=1;
FullTestLabels(num_test_pos+1:end)=-1;

for iPos = 1:num_test_pos
    fprintf('test pos: %d of %d.\n ',iPos,num_test_pos);
    load([cachePath '/' pos_file_list(test_p(iPos)).name]);
    FullTestData(iPos,:) = fVec(p1+1:p2);
end
for iNeg = 1:num_test_neg
    fprintf('test neg %d of %d.\n',iNeg,num_test_neg)
    load([cachePath '/' neg_file_list(test_n(iNeg)).name]);
    FullTestData(num_test_pos+iNeg,:)=fVec(p1+1:p2);
end


save([outPath  '/ABData_fold_' num2str(iFold) '.mat'],'FullData','FullLabels','FullTestData','FullTestLabels');


% Compute score of testing data, up to each feature; 
% AUC is also computed
for iFeature = 1:length(model)
[outScore]=adaboost('apply',FullTestData,model(1:iFeature));
auc=calc_auc(outScore,FullTestLabels');
aucTable(iFold,iFeature)=auc;
end
fileName = sprintf( '%s/Fold_%d_result.mat',outPath,iFold);
save(fileName,'aucTable','outScore','model');

% visulize active basis template;
vis_AB_template;
end

% Compute and plot AUC curve for active basis
meanAuc = mean(aucTable);
stdAuc = std(aucTable);
errorbar(1:nFeature,meanAuc,stdAuc,'r','linewidth',2)
axis([0 nFeature 0.5 1])

filename = sprintf('%s/AB_adaBoost_nFeature_%d_nTrain_%d_nTest_%d_nFold_%d',outPath,nFeature,num_train_pos+num_train_neg,num_test_pos+num_test_neg,nFold);
title(filename)
xlabel('number  of features');
ylabel('AUC');
saveas(gcf,[filename '.fig']);
save([filename '.mat'],'aucTable','outScore','model');

% Close the workers.
matlabpool close;
