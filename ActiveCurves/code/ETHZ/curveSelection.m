% maxI
% 1: orient
% 2: angle
% 3: length
% 4: x
% 5: y


% find the best response in shared M2^+ map
for iEl = 1:nLength*nAngle*norient
	[val ind]=max(reshape(pooledMap{iEl},sx*sy,1));
	pooledMax(iEl) = val;
	pooledMaxInd(iEl)=ind; 
end	
[val ind]=max(reshape(pooledMax,1,[]));

% find arg-max at arc level, and put into maxi array
[maxi(3),maxi(2),maxi(1)] = ind2sub(size(pooledMax),ind);
maxi = maxi-1;
ind = pooledMaxInd(maxi(3)+1,maxi(2)+1,maxi(1)+1);
[maxi(4),maxi(5)]=ind2sub([sx,sy],ind);
maxi(6)=val/(2*maxi(3)+1)/nimage;	
