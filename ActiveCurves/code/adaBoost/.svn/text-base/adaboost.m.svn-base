function [estimateclasstotal,model,dataweight]=adaboost(mode,datafeatures,dataclass_or_model, ...
                                            current_model,dataweight,itt)
% This function AdaBoost, consist of two parts a simpel weak classifier and
% a boosting part:
% The weak classifier tries to find the best treshold in one of the data
% dimensions to sepparate the data into two classes -1 and 1
% The boosting part calls the clasifier iteratively, after every classification
% step it changes the weights of miss-classified examples. This creates a
% cascade of "weak classifiers" which behaves like a "strong classifier"
%
%  Training mode:
%    [estimateclass,model]=adaboost('train',datafeatures,dataclass,itt)
%  Apply mode:
%    estimateclass=adaboost('apply',datafeatures,model)
% 
%  inputs/outputs:
%    datafeatures : An Array with size number_samples x number_features
%    dataclass : An array with the class off all examples, the class
%                 can be -1 or 1
%    itt : The number of training itterations
%    model : A struct with the cascade of weak-classifiers
%    estimateclass : The by the adaboost model classified data
%               
%  %% Example
%
%  example.m
%
%  Function is written by D.Kroon University of Twente (August 2010)

switch(mode)
    case 'train'
        % Train the adaboost model
        
        % Set the data class 
        dataclass=dataclass_or_model(:);
        if isempty(current_model)
            model=struct;
            t0=0;
        else
            model =current_model;
            t0 = length(current_model);
        end;
        
        % Weight of training samples, first every sample is even important
        % (same weight)
        if isempty(dataweight)
            dataweight=ones(length(dataclass),1)/length(dataclass);
        end
        % This variable will contain the results of the single weak
        % classifiers weight by their alpha
        estimateclasssum=zeros(size(dataclass));
        
        % Do all model training itterations
       for t=1:itt
            fprintf('Iteration: %d of %d\n',t+t0,itt+t0)
            % Find the best treshold to separate the data in two classes
            %[estimateclass,err,h] = WeightedThresholdClassifier(datafeatures,dataclass,D);
            WeightedThresholdClassifier;
            % Weak classifier influence on total result is based on the current
            % classification error
            alpha=1/2 * log((1-err)/max(err,eps));
            
            % Store the model parameters
            model(t+t0).alpha = alpha;
            model(t+t0).dimension=h.dimension;
            model(t+t0).threshold=h.threshold;
            model(t+t0).direction=h.direction;
            
            % We update D so that wrongly classified samples will have more weight
            dataweight =dataweight.* exp(-model(t+t0).alpha.*dataclass.*estimateclass);
            dataweight = dataweight./sum(dataweight);
            
            % Calculate the current error of the cascade of weak
            % classifiers
            estimateclasssum=estimateclasssum +estimateclass*model(t+t0).alpha;
            estimateclasstotal=sign(estimateclasssum);
            model(t+t0).error=sum(estimateclasstotal~=dataclass)/length(dataclass);
        %    if(model(t+t0).error==0), break; end
        end
    case 'apply' 
        % Apply Model on the test data
        
        % Get the Trained Adaboost Model
        model=dataclass_or_model;

        % Add all results of the single weak classifiers weighted by their alpha 
        estimateclasssum=zeros(size(datafeatures,1),1);
        for t=1:length(model);
            estimateclasssum=estimateclasssum+model(t).alpha*ApplyClassTreshold(model(t), datafeatures);
        end
        % If the total sum of all weak classifiers
        % is less than zero it is probablly class -1 otherwise class 1;
        % estimateclasstotal=sign(estimateclasssum);
        estimateclasstotal = estimateclasssum;
    otherwise
        error('adaboost:inputs','unknown mode');
end
end

%{
function [estimateclass,err,h] = WeightedThresholdClassifier(datafeatures,dataclass,dataweight)
% This is an example of an "Weak Classifier", it tries several tresholds
% in all data dimensions. The treshold which divides the data into two
% class with the smallest error is chosen as final treshold

% Set minimal treshold error in all dimensions to infinite
errdims = inf(1,size(datafeatures,2)); 
min_error = inf(1,1);

% Loop through the dimensions
nDim = size(datafeatures,2);
nSample = size(datafeatures,1);
whError = inf(nDim,1);
whThreshold = zeros(nDim,1,'single');
whDirection = zeros(nDim,1,'single');
parfor dim=1:nDim;
    % discard non-meaningful dimensions 
    if range(datafeatures(:,dim))<1e-6
        errdims(dim)=1e10; 
        continue;
    end
    if mod(dim,10000)==0
        fprintf('dimension %d of %d\n',dim,nDim)
    end
        [sorted_feature, order]= sort(datafeatures(:,dim),'descend');%direction ==1 
	sorted_class = dataclass(order);
        sorted_weight = dataweight(order);
        pos_weight = sorted_weight.*double((sorted_class~=1));
        neg_weight = sorted_weight.*double((sorted_class~=-1));
        test_error = zeros(1,nSample);
        test_error(1) = pos_weight(1)+sum(neg_weight(2:end));
        for iThresh = 2:nSample
           test_error(iThresh) = test_error(iThresh-1) + pos_weight(iThresh) - neg_weight(iThresh);
        end
        [val1 ind1]=min(test_error);
       	[val2 ind2] = max(test_error);
        if val1< 1-val2
        whError(dim) = val1;
	whDirection(dim)=1;
	whThreshold(dim) = mean([sorted_feature(ind),sorted_feature(min(ind+1,nSample))])+ 1e-10;
 
        else 
         whError(dim)=1-val2;
	 whDirection(dim) = 0; % meaning smaller value corresponds to +1(positive) output
         whThreshold(dim) = mean([sorted_feature(ind), sorted_feature(max(ind-1,0))])+1e-10;
         end
end
  [val ind]=min(whError); 
  h.direction = whDirection(ind);
  h.dimension = ind;
  h.threshold = whThreshold(ind);
  min_error = val;

  estimateclass=ApplyClassTreshold(h,datafeatures);
  err = min_error;
end
%}

function y = ApplyClassTreshold(h, x)
% Draw a line in one dimension (like horizontal or vertical)
% and classify everything below the line to one of the 2 classes
% and everything above the line to the other class.
if(h.direction == 1)
    y =  single(x(:,h.dimension) >= h.threshold);
else
    y =  single(x(:,h.dimension) < h.threshold);
end
y(y==0) = -1;
end
