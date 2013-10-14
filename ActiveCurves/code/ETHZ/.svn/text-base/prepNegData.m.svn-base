disp('preparing negative data');
load([inFolder '/Imgs/negImg.mat']);
Icolor = INeg;
nNeg  = numel(INeg);
for i = 1:nNeg
    INeg{i}=im2single(rgb2gray(Icolor{i}));
    Isize(i, :) = size(INeg{i});
end
sx = min(Isize(:, 1));
sy = min(Isize(:, 2));
Sx = sx*ones(1, nNeg);
Sy = sy*ones(1, nNeg);

M1map=cell(norient,1);
S2map=cell(nLength,nAngle,norient);
M2map=cell(nLength,nAngle,norient);
for (o = 1 : norient)
    M1map{o} = -1e10*ones(sx, sy,'single');
end
for iEL = 1:norient*nAngle*nLength
    S2map{iEL} = NEGMAX*ones(sx, sy,'single');
    M2map{iEL} = NEGMAX*ones(sx, sy,'single');
end



for i = 1:nNeg
    S1map = applyfilterfftsame(INeg(i), allfilter);  % filter training images
    ClocalNormalizeSingle(Sx(i),Sy(i),norient,h,localHalfx,localHalfy,...
        S1map,thresholdFactor);
    CsigmoidSingle(1,Sx,Sy,norient,saturation,S1map);
    outFileName = [workFolder  '/S1_neg_' num2str(i,'%04d') 'ori' '.mat'];
    save(outFileName,'S1map');
    
    CMax1Single(1,norient, S1map, M1map, Lrange, Orange, Sx, Sy);
    
    CSum2Single(1, norient, nLength,(nAngle-1)/2, M1map, S2map, sx, sy, h,lambda,logZ);
   
    
    outFileName = [workFolder  '/S2_neg_' num2str(i,'%04d') 'ori' '.mat'];
    save(outFileName,'S2map');

end
