%% assume bdImage is the image bdImage;
%%bdCell
IDtemp = zeros(size(I));

nTest = length(bbImgCell);

Wx = ones(1,1)*wx;
Wy = ones(1,1)*wy; 

imaxi = zeros(1, 6);
for iTest = 1:nTest
   dsym = zeros(wx,wy,'single');
    S1map = applyfilterfftsame(bbImgCell(iTest), allfilter);  % filter training images
    ClocalNormalizeSingle(wx,wy,norient,h,localHalfx,localHalfy,...
        S1map,thresholdFactor);
    CsigmoidSingle(1,Wx,Wy,norient,sat,S1map);
    M1map=cell(norient,1);
    for o = 1 : norient
        M1map{o} = -1e10*ones(wx, wy,'single');
    end

    S2map = cell(nLength,nAngle,norient);
    M2map = cell(nLength,nAngle,norient);
    for iAng = 1:nAngle
        for iOri = 1:norient
            for iLen = 1:nLength
            S2map{iLen,iAng,iOri} = -1e10*ones(wx, wy,'single');
            M2map{iLen,iAng,iOri} = -1e10*ones(wx, wy,'single');
            end
        end
    end
    CMax1Single(1,norient, S1map, M1map, Lrange, Orange, Wx, Wy);
    
    CSum2Single(1, norient, nLength,(nAngle-1)/2, M1map, S2map, wx, wy, h,lambda,logZ);
    
    CMax2Single(1, norient, nLength, nAngle, S2map, M2map, Lrange2, Orange2, Wx, Wy);
    for iCurve = 1:size(selCurves,1);
        maxi = selCurves(iCurve,:);
    CBackMax2Single(1, norient, nLength, nAngle, S2map, h, Lrange2, Orange2, 1, Wx, Wy, maxi, imaxi);
    CdrawCurveSingle(1, norient, nLength, nAngle, S1map, wx, wy, h, Lrange, Orange, maxi, imaxi,  dsym, allsymbol(1,:));
    end
    %% paste the template back;
    bb = BBs(iTest,5:8);
    bb = [max(1,bb(1)),max(1,bb(2)),min(size(I,1),bb(3)),min(size(I,2),bb(4))];
    bb = round(bb);
    if(bb(3)-bb(1))<0 ||(bb(4)-bb(2))<0
        continue;
    end
    dsym = imresize(dsym,[bb(3)-bb(1)+1, bb(4)-bb(2)+1]);
    I(bb(1):bb(3),bb(2):bb(4))=max(I(bb(1):bb(3),bb(2):bb(4)),dsym);
    IDtemp(bb(1):bb(3),bb(2):bb(4))=dsym;
    
%     
%     bb =bbs(itest);
%     bd = [bb.top/bb.scale/bb.asp, bb.left/bb.scale, wx/bb.scale/bb.asp, wy/bb.scale];
%     bd = round(bd);
%     bd_img=[1 1 size(i,1) size(i,2)];
%     bd_int = [max(bd_img(1:2),bd(1:2)) min(bd_img(1:2)+bd_img(3:4)-1,bd(1:2)+bd(3:4)-1)];
%     bd_int(3:4)=bd_int(3:4)-bd_int(1:2)+1;
%     
%     dsym = imresize(dsym,[bd(3), bd(4)]);
%     sym_bd = [bd_int(1:2)-bd(1:2)+1 bd_int(3:4)];
%     i(bd_int(1):bd_int(1)+bd_int(3)-1,bd_int(2):bd_int(2)+bd_int(4)-1)=max(...
%         i(bd_int(1):bd_int(1)+bd_int(3)-1,bd_int(2):bd_int(2)+bd_int(4)-1),...
%         dsym(sym_bd(1):sym_bd(1)+sym_bd(3)-1,sym_bd(2):sym_bd(2)+sym_bd(4)-1));
%     
%     idtemp(bd_int(1):bd_int(1)+bd_int(3)-1,bd_int(2):bd_int(2)+bd_int(4)-1) = ...
%         dsym(sym_bd(1):sym_bd(1)+sym_bd(3)-1,sym_bd(2):sym_bd(2)+sym_bd(4)-1);
end
