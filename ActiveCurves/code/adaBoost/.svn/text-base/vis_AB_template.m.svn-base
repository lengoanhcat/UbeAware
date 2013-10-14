%% visualize active basis selected by adaboost
n_sel_features = length(model);



sel_feature_index = zeros(n_sel_features,6);
% convert index to geometric attributes
for i_feature = 1:n_sel_features
    dim = model(i_feature).dimension;
    y = ceil(dim/sx/nOrient);
    residual_dim = dim-(y-1)*sx*nOrient;  
    i_orient = ceil(residual_dim/sx);
    x= residual_dim - (i_orient-1)*sx;
    sel_feature_index(i_feature,:) = [ i_orient-1,0,0,x,y,model(i_feature).alpha];
end 


% draw the template, 
S1maps = cell(nOrient,1);
[allfilter, allsymbol] = makefilter(scale, nOrient);  % generate Gabor filters
h = (size(allfilter{1}, 1)-1)/2; % half size of Gabor
C = corr(allfilter, epsilon);    % generate the inhibition maps

sym = zeros(sx,sy,'single');
Asym = cell(1);
Asym{1}=zeros(sx,sy,'single');

for i_feature = 1:n_sel_features
	for i_orient = 1:nOrient
		S1map{i_orient}=model(i_feature).alpha*ones(sx,sy,'single'); 
	end
	maxi = sel_feature_index(i_feature,:);
	imaxi = maxi;
	ds_respVec = zeros(2*imaxi(3)+1,1,'single');
	CBackM1andInhibition(1, nOrient,1,1,S1map, C,...
	                sx,sy, h, Lrange, Orange,...
	                maxi, imaxi, ds_respVec, sym, Asym, allsymbol(1, :),0);
	img = sym/max(sym(:));
	imgName= sprintf('%s/AB_template_Fold_%d_feature_%d.png',outPath,iFold,i_feature);
	imwrite(1-img,imgName);
end
