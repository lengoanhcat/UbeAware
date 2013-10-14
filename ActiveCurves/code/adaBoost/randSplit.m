function [train_ind test_ind]=randSplit(numImage,trainPortion);
index = randperm(numImage);
nTrain = floor(numImage*trainPortion);
train_ind = index(1:nTrain);
test_ind = index(nTrain+1:end);