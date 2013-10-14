function BBs= testSingleImage(I,Iparam,DParam,HParam)
IColor = I;
if(size(I,3)==3)
    I = rgb2gray(I);
end
disp(['begin testing image: ' Iparam.imgName] )

%% Filter Param
allfilter = HParam.allfitler;
norient = HParam.norient;
h = HParam.h;
localHalfx= HParam.localHalfx;
localHalfy=HParam.localHalfy;
Lrange = HParam.Lrange;
Orange = HParam.Orange;
thresholdFactor = HParam.thresholdFactor;
sat = HParam.sat;
allsymbol = HParam.allsymbol;

%% Detector Parameter
selCurves=DParam.selCurves;

buff_path = DParam.buffer_path;
nLength = DParam.nLength;
nAngle = DParam.nAngle;
lambda = DParam.lambda;
logZ = DParam.logZ;
Lrange2 = DParam.Lrange2;
Orange2 = DParam.Orange2;
wx = DParam.wx;
wy = DParam.wy;
maxWindowNum = DParam.maxWindowNum;
pad = DParam.image_pad;
LRWeight = DParam.LRWeight;



sx = size(I,1); sy = size(I,2);
% decide minimum and maximum image scale
maxS = wy*6/sy;
minS = 5/6*wy/sy;
dS = 1.1;
minAsp = 0.8;maxAsp = 1.2;
dAsp = 0.1;
n = round(log(maxS/minS)/log(dS)); % number of scales
s = minS*(dS.^(0:n)); % scales
asp = minAsp:dAsp:maxAsp; % aspect ratios

nCurve = size(selCurves,1); % number of curves used in template
nDeform = (Orange2*2+1); % number of deformations used


cateInd = strfind(Iparam.path,'/');
cacheFile = [DParam.outPath '/' Iparam.path(cateInd(end)+1:end) '_' Iparam.imgName(1:end-3) '_scoremap.mat'];

if(~ exist(cacheFile,'file'))
    scoreMaps = cell(length(s),length(asp));
    for is = 1:length(s)
        disp(['scale ' num2str(is) ' of ' num2str(length(s))])
        for iAsp = 1:length(asp)
            sx = size(I,1); sy = size(I,2);
            if( round(sy*s(is)+2*h)<=wy || round(sx*s(is)*asp(iAsp)+2*h)<=wx)
                % skip this deformation
                scoreMaps{is,iAsp} = -1e10*ones(1,1);
                continue;
            end

            img = imresize(I,round([sx*s(is)*asp(iAsp) sy*s(is)]));
            Iscale{1}=padarray(img,[pad pad],'symmetric');
            sx = size(Iscale{1},1); sy = size(Iscale{1},2);
            S1map = applyfilterfftsame(Iscale, allfilter);  % filter training images
            ClocalNormalizeSingle(sx,sy,norient,h,localHalfx,localHalfy,...
                S1map,thresholdFactor);
            CsigmoidSingle(1,sx,sy,norient,sat,S1map);

            M1map=cell(norient,1);
            for o = 1 : norient
                M1map{o} = -1e10*ones(sx, sy,'single');
            end
            CMax1Single(1,norient, S1map, M1map, Lrange, Orange, sx, sy);

            %% S2 maps and M2 maps
            S2map=cell(nCurve,nDeform);
            M2map=cell(nCurve,1);
            for iCurve = 1:nCurve
                M2map{iCurve,1}=-1e10*ones(sx, sy,'single');
                for iDeform = 1:nDeform
                    S2map{iCurve,iDeform}=-1e10*ones(sx, sy,'single');
                end
            end
            %tic
            CTestS2M2Single(norient, nLength,nAngle,Lrange2,Orange2,...
			    nCurve,selCurves,M1map, S2map,M2map, ...
                            sx, sy, h,single(LRWeight));
            %disp(['S2M2map: '  num2str(toc)])
            %% begin detection by sliding windows
            s_map = slideWindow(M2map,selCurves,wx,wy);
            scoreMaps{is,iAsp}=s_map;
        end
    end

    save(cacheFile,'scoreMaps','I','DParam','HParam');
else
    load(cacheFile);
end
%% process the score maps to report bd
bbImgCell= cell(1,1);
sx = size(I,1); sy = size(I,2);
iBB =0 ;
while(1)
    iBB = iBB + 1;
    % select max
    max_scores =zeros(size(scoreMaps));
    max_inds = zeros(size(scoreMaps));
    for is = 1:length(s)
        for iAsp = 1:length(asp)
            [max_scores(is,iAsp) max_inds(is,iAsp)] = max(scoreMaps{is,iAsp}(:));
        end
    end
    [max_val max_ind]=max(max_scores(:));

    % fill the attributes of BB structure
    bb.score = max_val;
    [is iAsp]=ind2sub(size(scoreMaps),max_ind);
    bb.scale = s(is);
    bb.wx = wx;
    bb.wy = wy;
    bb.asp = asp(iAsp);
    [bb.top bb.left] = ind2sub(size(scoreMaps{is,iAsp}),max_inds(is,iAsp));
    bb.top = bb.top-pad;
    bb.left = bb.left-pad;
    BBs(iBB) = bb;

    % crop bbimage 
    bd = [bb.top/bb.scale/bb.asp, bb.left/bb.scale, (bb.top+wx)/bb.scale/bb.asp, (bb.left+wy)/bb.scale];
    bd = round(bd);
    window = subarraySingle(I, bd(1), bd(3), bd(2), bd(4), 1);
    window = imresize(window,[wx wy]);
    bbImgCell{iBB}=window;

    if(bb.score+LRWeight(end)<0 || iBB >=maxWindowNum)
        break;
    end
    % small bouding box in original image coordinate
    % The padding on training image is 10
    sTop = (bb.top+10)/bb.scale/bb.asp;
    sBot = (bb.top +10+ wx-2*10)/bb.scale/bb.asp;
    sLeft = (bb.left+10)/bb.scale;
    sRight = (bb.left+10+wy-2*10)/bb.scale;
    % inhibition
    for is = 1:length(s)
        for iAsp = 1:length(asp)
	    % smalBB size in original image coord
            [height width]=size(scoreMaps{is,iAsp});
            swx = (wx-20)/s(is)/asp(iAsp);
            swy = (wy-20)/s(is);
	    
            %% inhibite small BB in original image coord
            inh_top = sTop - swx; inh_left = sLeft - swy;
            inh_bot = sBot ; inh_right = sRight ;
	     
            % to target image coord, 
	    % 10 is the padding on object BB
            % pad is the padding on image boundary
            ratio = s(is)*asp(iAsp);
            row0 = round(max(1,(inh_top)*ratio-10+pad));
            row1 =round( min(height,(inh_bot)*ratio-10+pad));
            ratio = s(is);
            col0=round(max(1,(inh_left)*ratio-10+pad));
            col1 =round( min((inh_right)*ratio-10+pad,width));
            scoreMaps{is,iAsp}(row0:row1,col0:col1) = -1e10;
        end
    end
end

[res score]=evalBBs(BBs,Iparam,10);
%% draw individual template
drawBBTemplate;
imwrite(IDtemp,[DParam.outPath '/' Iparam.path(cateInd(end)+1:end) '_' Iparam.imgName(1:end-3) '_templates.png']);
imagesc(I);


%% draw bounding box

colormap(gray);
hold on;
bb=[];
% ground truth
if isfield(Iparam,'anno')
    for iGT = 1:size(Iparam.anno,1)
        bb= Iparam.anno(iGT,:);
        plot(bb([1 3 3 1 1]),bb([2 2 4 4 2]),'b-','linewidth',2);
    end
end
% predicted;
for iBB = 1:length(BBs)
    bb(1) = (BBs(iBB).left+10)/ BBs(iBB).scale;
    bb(2) = (BBs(iBB).top+10)/ BBs(iBB).scale/BBs(iBB).asp;
    bb(3) = (BBs(iBB).wy-20)/BBs(iBB).scale+bb(1);
    bb(4) = (BBs(iBB).wx-20)/BBs(iBB).scale/BBs(iBB).asp+bb(2);
    if( res(iBB)==1)
        plot(bb([1 3 3 1 1]),bb([2 2 4 4 2]),'g-','linewidth',2);
        text(bb(1)+2,bb(2)+8,num2str(score(iBB)),'Color','g')
        text(bb(1)+2,bb(4)-15,num2str(BBs(iBB).asp),'Color','g')
    else
        %plot(bb([1 3 3 1 1]),bb([2 2 4 4 2]),'r-','linewidth',2);
        %text(bb(1)+2,bb(2)+8,num2str(score(iBB)),'Color','r')
        %text(bb(1)+2,bb(4)-15,num2str(BBs(iBB).asp),'Color','r')
    end
end
hold off;
axis image;

saveas(gcf,[DParam.outPath '/' Iparam.path(cateInd(end)+1:end) '_' Iparam.imgName(1:end-3) '_result.png']);
close ;
end
