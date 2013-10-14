function BBs= testSingleImageV2(I,Iparam,DParam,HParam)
% compared with V1, the main change is window reporting strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%minS = 1/8; maxS = 1/2;
maxS = wy*6/sy;
minS = 5/6*wy/sy;
dS = 1.1;
minAsp = 0.8;maxAsp = 1.2;
dAsp = 0.1;
n = round(log(maxS/minS)/log(dS));
s = minS*(dS.^(0:n));
asp = minAsp:dAsp:maxAsp;


nCurve = size(selCurves,1);
nDeform = (Orange2*2+1);


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
            %% begin detection
            s_map = slideWindow(M2map,selCurves,wx,wy);
            scoreMaps{is,iAsp}=s_map;
        end
    end

    save(cacheFile,'scoreMaps','I','DParam','HParam');
else
    load(cacheFile);
end
%% process the score maps to report bd

sx = size(I,1); sy = size(I,2);
threshold = -LRWeight(end);
obj_top=[];obj_bot=[];obj_left=[];obj_right=[];
crop_top=[];crop_bot=[];crop_left=[];crop_right=[];
obj_bb=[];
crop_bb=[];
score = [];
bb_scale =[];
bb_asp = [];
% find all possible boxes that is above threshold;
for is = 1:length(s)
    for iAsp = 1:length(asp)
        ind = find(scoreMaps{is,iAsp}>threshold);
        [row, col]=ind2sub(size(scoreMaps{is,iAsp}),ind);
        obj_bb = [ obj_bb;...
                  (row-pad +10)/s(is)/asp(iAsp),(col-pad +10)/s(is),...  
                  (row-pad +wx -10)/s(is)/asp(iAsp),(col-pad +wy -10)/s(is)];
              
        crop_bb = [ crop_bb;...
                  (row-pad )/s(is)/asp(iAsp),(col-pad)/s(is),...  
                  (row-pad +wx)/s(is)/asp(iAsp),(col-pad +wy)/s(is)];    
        score = [score;   scoreMaps{is,iAsp}(ind)]; 
        bb_scale = [bb_scale; s(is)*ones(length(ind),1)];
        bb_asp = [bb_asp; asp(iAsp)*ones(length(ind),1)];
    end
    pic =nms([obj_bb,score],0.5);
    obj_bb=obj_bb(pic,:);
    crop_bb=crop_bb(pic,:);
    score = score(pic,:);
    bb_scale = bb_scale(pic,:);
    bb_asp = bb_asp(pic,:);
end
% non-maximum supression
pick =nms([obj_bb,score],0.5);
all_bb =[round(obj_bb), round(crop_bb),bb_scale,bb_asp,score];
all_bb = all_bb(pick,:);
if(size(all_bb,1)>5)
    all_bb=all_bb(1:5,:);
end
% crop images and draw templates
n_bb = size(all_bb,1);
bbImgCell= cell(n_bb,1);
for i = 1:n_bb
    window = subarraySingle(I, all_bb(i,1+4), all_bb(i,3+4), all_bb(i,2+4), all_bb(i,4+4), 1);  
    window = imresize(window,[wx wy]);
    bbImgCell{i}=window;
end

BBs = all_bb;

%drawBBTemplateV2;
%imwrite(IDtemp,[DParam.outPath '/' Iparam.path(cateInd(end)+1:end) '_' Iparam.imgName(1:end-3) '_templates.png']);

imagesc(I);


%% draw bounding box

[res score]=evalBBsV2(BBs,Iparam);
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
for iBB = 1:n_bb
    bb = double(BBs(iBB,[2 1 4 3]));
    if( res(iBB)==1)
        plot(bb([1 3 3 1 1]),bb([2 2 4 4 2]),'g-','linewidth',2);
        text(bb(1)+2,bb(2)+8,num2str(BBs(iBB,end)),'Color','g')
        text(bb(1)+2,bb(4)-15,num2str(BBs(iBB,end-1)),'Color','g')
    else
        plot(bb([1 3 3 1 1]),bb([2 2 4 4 2]),'r-','linewidth',2);
        text(bb(1)+2,bb(2)+8,num2str(BBs(iBB,end)),'Color','r')
        text(bb(1)+2,bb(4)-15,num2str(BBs(iBB,end-1)),'Color','r')
    end
end
hold off;
axis image;
saveas(gcf,[DParam.outPath '/' Iparam.path(cateInd(end)+1:end) '_' Iparam.imgName(1:end-3) '_result.png']);
close ;
end