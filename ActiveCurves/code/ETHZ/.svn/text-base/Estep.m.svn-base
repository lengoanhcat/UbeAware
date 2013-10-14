for iImg = 1:nimage
disp([ 'Estep to ' num2str(iImg) '/' num2str(nimage)]);
inFileName = [workFolder  '/S1_'  num2str(iImg,'%04d') 'ori' '.mat'];
load(inFileName);

disp(['get M1 maps']);
CMax1Single(1,norient, S1map, M1map, Lrange, Orange, Sx, Sy);

disp(['get S2 maps of line segments']);
CSum2Single(1, norient, nLength,(nAngle-1)/2, M1map, S2map, sx, sy, h,lambda,logZ);

disp(['get M2 score']);
CMax2Single(1, norient, nLength, nAngle, S2map, M2map, Lrange2, Orange2, Sx, Sy);

for iCluster = 1:numCluster
	selCurves = clusterTemplate{iCluster};
	score = 0;
	for iCurve = 1:size(selCurves,1)
		maxi = selCurves(iCurve,:);
		score = score + M2map{maxi(3)+1,maxi(2)+1,maxi(1)+1}(maxi(4),maxi(5));
	end
	scoreMat(iImg,iCluster)=score;
end

disp(num2str(scoreMat(iImg,:)));
[val ind]=max(scoreMat(iImg,:));
disp(['Image ' num2str(iImg)  ' is assigned to cluster ' num2str(ind)])
end
[val cluster] = max(scoreMat');
