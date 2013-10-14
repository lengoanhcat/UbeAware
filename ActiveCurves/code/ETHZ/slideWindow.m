function score_map = slideWindow(M2map,selCurves,wx,wy)
% wx: window size , x
% wy: window size, y
    [sx sy] = size(M2map{1});
    score_map = zeros(sx-wx,sy-wy,'single');
    for iCurve = 1:size(selCurves,1)
        maxi = selCurves(iCurve,:);
        m2_map = M2map{iCurve};
        row0=max(1,maxi(4));
        row1=min(sx-wx+maxi(4)-1,sx);
        col0=max(1,maxi(5));
        col1=min(sy-wy+maxi(5)-1,sy);
        score_map = score_map + m2_map(row0:row1,col0:col1);
    end
end
