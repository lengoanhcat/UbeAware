function cropData(inPath,outPath,aligned_height,marg,prop,PNRatio)
%
% prop: proportion of data used for training
global FaceDB;

outPath = [outPath '/Imgs']
% clean the output folder
if(~exist(outPath))
	mkdir(outPath)
else
	delete([outPath '/*.*'])
end

flg_exist_bb = true; % default bounding boxes data are included

maxWidth = 0; % pre-allocated for selecting the image width
if ~isempty(dir([inPath '/*.groundtruth']))
    inFiles = dir([inPath '/*.groundtruth']);
elseif ~isempty(dir([inPath '/*.bbx']))
    inFiles = dir([inPath '/*.bbx']);
elseif ~isempty(dir([inPath '/*.jpg']))
    inFiles = dir([inPath '/*.jpg']);
    flg_exist_bb = false;
end

anno=[];

numImg = length(inFiles);
indTrain = 1:round(prop*numImg);
indTest = round(prop*numImg)+1:numImg;


%% begin cropping positive training data
totalInst=1;
aspHist = [];
% read bounding box and create aspect ratio statistics
for i = 1:length(indTrain)
    inFiles(indTrain(i)).name;
    data = load([ inPath '/' inFiles(indTrain(i)).name]);
    nInst = size(data,1);
    if FaceDB ~= 1
        pos = strfind(inFiles(indTrain(i)).name,'_');        
    else
        pos = strfind(inFiles(indTrain(i)).name,'.');        
    end    
    pos = pos(end);        
    
    imgName = [ inPath '/' inFiles(indTrain(i)).name(1:pos-1) '.jpg'];    
    if exist(imgName,'file') ~= 2
        imgName = [ inPath '/' inFiles(indTrain(i)).name(1:pos-1) '.JPG'];
    end

    for iInst=1:nInst
        anno(totalInst).bd = data(iInst,:);
        anno(totalInst).imgName = imgName;
        %% crop object instance
        width = data(iInst,3)-data(iInst,1);
        height = data(iInst,4)-data(iInst,2);
        aspHist(totalInst)=width/height;
        anno(totalInst).width = width;
        anno(totalInst).height = height;
        anno(totalInst).instName = [ outPath '/' inFiles(indTrain(i)).name(1:pos-1) '_inst' num2str(iInst) '.bmp'];
        maxWidth = max(maxWidth,aligned_height/height*width);
        totalInst = totalInst+1;        
    end
end
% decide bounding box size
meanAsp = geomean(aspHist);
aligned_width = aligned_height*meanAsp;
IPos = cell(1);
% begin cropping image
for iInst = 1:length(anno)
    bd = anno(iInst).bd;
    bd = round(bd);
    img = imread(anno(iInst).imgName);
    if size(img,3)==1
    	img = gray2rgb(img);
    end
    scale_h = aligned_height/anno(iInst).height;
    pady = round(marg /scale_h);
    scale_w = aligned_width/anno(iInst).width;
    padx = round(marg/scale_w);
    window = subarray(img, bd(2)-pady, bd(4)+pady, bd(1)-padx, bd(3)+padx, 1);
    window = imresize(window,[aligned_height+2*marg aligned_width+2*marg]);
    IPos{iInst}=window;
    imwrite(window,anno(iInst).instName);
end
save([outPath '/posImg.mat'],'IPos');
%% crop negative data
numNeg = round(1/PNRatio);
% For each cropped boundingbox, we randomly create numNeg bounding bodx
% as natrual image
INeg = cell(1);
totalNeg = 0;
for iInst = 1:length(anno)
    bd = anno(iInst).bd;
    bd=round(bd);
    img = imread(anno(iInst).imgName);
    if size(img,3)==1
    	img = gray2rgb(img);
    end
    [h w d]=size(img);
    for iNeg = 1:numNeg
       validBB=false;
       while ~validBB      
       negBD = 0.1+0.8*rand(4,1);
       negBD(1)=negBD(1)*w;
       negBD(2)=negBD(2)*h;
       negBD(3) = negBD(3).*(w-negBD(1));
       negBD(4) = negBD(3)/meanAsp;
       negBD(4) = min(negBD(4),h-negBD(2)-1);
       negBD = round(negBD);
       negBD(3) = negBD(3)+negBD(1)-1;
       negBD(4) = negBD(2)+negBD(4)-1;
       if(negBD(3)<w && negBD(4)<h)
          validBB=true;
       end
       end
        totalNeg = totalNeg + 1;
        window = img(negBD(2):negBD(4),...
                     negBD(1):negBD(3),:);
        window = imresize(window,[aligned_height+2*marg aligned_width+2*marg]);
        INeg{totalNeg}=window;
        imwrite(window,...
               [anno(iInst).instName(1:end-4) 'neg' num2str(iNeg) '.bmp']);
    end
end
save([outPath '/negImg.mat'],'INeg');
%% save testing data index
testImg = [];
for i = 1:length(indTest)
    data = load([ inPath '/' inFiles(indTest(i)).name]);
    testImg(i).nInst = size(data,1);
    if FaceDB == 1
        pos = strfind(inFiles(indTest(i)).name,'.');
    else 
        pos = strfind(inFiles(indTest(i)).name,'_');
    end
    pos = pos(end);
    testImg(i).path = inPath;
    testImg(i).imgName = [ inFiles(indTest(i)).name(1:pos-1) '.jpg'];
	anno = [];
    for iInst=1:nInst
        anno(iInst).bd = data(iInst,:);
    end
	testImg(i).anno = data;
end
save([outPath '/testIndex.mat'],'testImg')


