
% This is an example of an "Weak Classifier", it tries several tresholds
% in all data dimensions. The treshold which divides the data into two
% class with the smallest error is chosen as final treshold

% Set minimal treshold error in all dimensions to infinite
%errdims = inf(1,size(datafeatures,2));
min_error = inf(1,1);

% Loop through the dimensions
nDim = size(datafeatures,2);
nSample = size(datafeatures,1);
whError = inf(nDim,1);
whThreshold = zeros(nDim,1,'single');
whDirection = zeros(nDim,1,'single');
step = 4e5;

for batch_num = 1:step:nDim;
	fprintf('to %d of %d\n',batch_num, nDim);
	batch_end = min(batch_num+step,nDim);
	crop_nDim = batch_end-batch_num+1;
	crop_datafeatures = datafeatures(:,batch_num:batch_end);
        crop_whError = inf(crop_nDim,1);
        crop_whThreshold = zeros(crop_nDim,1,'single');
	crop_whDirection = zeros(crop_nDim,1,'single');

parfor dim=1:crop_nDim;
    % discard non-meaningful dimensions 
    if range(crop_datafeatures(:,dim))<1e-6
        %errdims(batch_num-1+dim)=1e10;
        continue;
    end
        [sorted_feature, order]= sort(crop_datafeatures(:,dim),'descend');%direction ==1 
        sorted_class = dataclass(order);
        sorted_weight = dataweight(order);
        pos_weight = sorted_weight.*double((sorted_class~=1));
        neg_weight = sorted_weight.*double((sorted_class~=-1));
        test_error = zeros(1,nSample);
        %{ 
	test_error(1) = pos_weight(1)+sum(neg_weight(2:end));
        for iThresh = 2:nSample
           test_error(iThresh) = test_error(iThresh-1) + pos_weight(iThresh) - neg_weight(iThresh);
        end
	verify_test_error = test_error;
        %}

	acc_pos_weight = cumsum(pos_weight);
	acc_reverse_neg_weight = cumsum(neg_weight(nSample:-1:1));
	test_error(1:nSample-1)= acc_pos_weight(1:nSample-1)+acc_reverse_neg_weight(nSample-1:-1:1);
        test_error(nSample)=acc_pos_weight(nSample);

	%{
	diff=sum(abs(verify_test_error-test_error));
	fprintf('error: %f',diff);
        %}

        [val1 ind1]=min(test_error);
        [val2 ind2] = max(test_error);
        if val1< 1-val2
        crop_whError(dim) = val1;
        crop_whDirection(dim)=1;
        crop_whThreshold(dim) = mean([sorted_feature(ind1),sorted_feature(min(ind1+1,nSample))])+ 1e-10;
        else
        crop_whError(dim)=1-val2;
        crop_whDirection(dim) = 0; % meaning smaller value corresponds to +1(positive) output
        crop_whThreshold(dim) = mean([sorted_feature(ind2), sorted_feature(max(ind2-1,1))])+1e-10;
         end
end % parfor
	whError(batch_num:batch_end) = crop_whError;
	whDirection(batch_num:batch_end)=crop_whDirection;
	whThreshold(batch_num:batch_end)=crop_whThreshold;
end

  [val ind]=min(whError);
  h.direction = whDirection(ind);
  h.dimension = ind;
  h.threshold = whThreshold(ind);
  min_error = val;

  %estimateclass=ApplyClassTreshold(h,datafeatures);
  if(h.direction == 1)
    y =  single(datafeatures(:,h.dimension) >= h.threshold);
  else
    y =  single(datafeatures(:,h.dimension) < h.threshold);
  end
  y(y==0) = -1;
  estimateclass = y;
  err = min_error;
