function testCategoryV2(name,buff_path,run_id)
addpath(genpath([pwd '/LargeProcessing']));
mex CTestS2M2Single.c
mex CMax1Single.c
mex CSum2Single.c
mex CMax2Single.c
mex CdrawCurveSingle.c
mex ClocalNormalizeSingle.c
mex CsigmoidSingle.c;
mex CBackMax2Single.c

% test images over one category
%% load parameter for that category
global FaceDB;
cd(['../../data/cache/ETHZ/extractData' num2str(run_id) '/' name]);
param;
cd ('./Imgs')
load testIndex;
cd ('../../../../../../code/ETHZ');
modelName = ['../../results/ETHZ/model' num2str(run_id) '/' name '_template.mat'];
load(modelName);

[allfilter, allsymbol] = makefilter(scale, norient);  % generate Gabor filters
h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor

% fillout DParam and HParam
HParam.allfitler=allfilter ;
HParam.allsymbol = allsymbol;
HParam.norient=norient;
HParam.h=h;
HParam.localHalfx=localHalfx;
HParam.localHalfy=localHalfy;
HParam.Lrange=Lrange;
HParam.Orange=Orange;
HParam.thresholdFactor = thresholdFactor;
HParam.sat=6;





DParam.selCurves=selCurves; % selected arcs
DParam.buffer_path=buff_path; % place for storing bounding-boxed images
DParam.nLength=nLength; % length of curves
DParam.nAngle=nAngle; % # of curvatures
DParam.lambda=lambda; % the coefficient
DParam.logZ=logZ ; 
DParam.Lrange2=Lrange2; % location deformation range
DParam.Orange2=Orange2; % orintation deformation range
DParam.wx=sx; % template size
DParam.wy=sy; % template size
DParam.maxWindowNum=5; % number of reported windows
DParam.drawBD= true; % if draw the boundingbox
DParam.drawTmplate = false; % if draw the deformed template
DParam.outPath = ['../../results/ETHZ/result' num2str(run_id) '/' name]; % result output path
if FaceDB == 1
    DParam.dataPath = ['../../data/ETHZ/Swans/']; % testing data path, should be the directory of ETHZ dataset
else 
    DParam.dataPath = ['../../data/ETHZ/']; % testing data path, should be the directory of ETHZ dataset
end
DParam.image_pad = 20; % image padding for detecing object near the image boundary.
DParam.LRWeight = LRWeight; % weight of each feature
if ~exist(DParam.outPath,'dir');
    mkdir(DParam.outPath);
end

%% load testing image.

% positive images
num_pos_inst = 0;
posRes =[];
posScore=[];
for iFile = 1:length(testImg)
    num_pos_inst = num_pos_inst + size(testImg(iFile).anno,1);
    img = im2single(imread([testImg(iFile).path '/' testImg(iFile).imgName]));
    BBs=testSingleImageV2(img,testImg(iFile),DParam,HParam);
    testImg(iFile).BBs = BBs;
    [res score]=evalBBsV2(BBs,testImg(iFile));
    posRes = [posRes; res];
    posScore = [posScore;score];
end
S.tot_pos = length(testImg);
% negtive images
negRes =[];
negScore=[];
folders = dir(DParam.dataPath);
for iFolder = 1:length(folders)
    if( strcmp(folders(iFolder).name,'.') ||...
            strcmp(folders(iFolder).name,'..') ||...
            strcmp(folders(iFolder).name,'horse')||...
            strcmp(folders(iFolder).name,name))
        continue;
    end
    imgs = dir([DParam.dataPath '/' folders(iFolder).name '/*.jpg']);
    for iImg = 1:length(imgs)
        iFile = iFile + 1;
        testImg(iFile).path = [DParam.dataPath '/' folders(iFolder).name];
        testImg(iFile).imgName = imgs(iImg).name;
        img = im2single(imread([testImg(iFile).path '/' testImg(iFile).imgName]));
        BBs=testSingleImageV2(img,testImg(iFile),DParam,HParam);
        testImg(iFile).BBs=BBs;
        [res score]=evalBBsV2(BBs,testImg(iFile));
        negRes = [negRes;res];
        negScore = [negScore;score];
    end
end
close all;


S.tot_test_imgs = iFile;
S.tot_corr = num_pos_inst;
S.tot_neg = S.tot_test_imgs-S.tot_pos;
S.corr_scores = [posScore(posRes ==1)]';

S.false_scores = [posScore(posRes ==0);negScore]';

save(['../../results/ETHZ/result' num2str(run_id) '/' name '_SV2.mat'],'S','posRes','posScore','negRes','negScore');
OutputStatistics(S,['../../results/ETHZ/' num2str(run_id) '/' name ]);
close;