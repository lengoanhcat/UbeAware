pcode curveSelection;

addpath liblinear/matlab

disp(['get constants']);


disp(['get C-codes']);
mex CsigmoidSingle.c
mex ClocalNormalizeSingle.c
mex CMax1Single.c;  % compute the MAX1 map
mex CSum2Single.c;  % compute the SUM2 map
mex CMax2Single.c;  % compute the MAX2 score
mex CBackMax2Single.c;
mex CBackM1andInhibition.c % inhibit the SUM1 map
dataFolder = ['../../data/cache/ETHZ/extractData' num2str(run_id)];
categories = dir(dataFolder);


workFolder =[ './working' num2str(run_id)];
mkdir(workFolder);


outFolder = ['../../results/ETHZ/model' num2str(run_id)];
if ~exist(outFolder)
   mkdir(outFolder)
end
mkdir(outFolder)
for indCate = 1:length(categories)
    if( strcmp(categories(indCate).name ,'.') ...
            || strcmp(categories(indCate).name , '..')...
            || strcmp(categories(indCate).name , '.DS_Store'))
        continue;
    end
    disp('get parameter file')
	inFolder = [ dataFolder '/' categories(indCate).name ] ;  % image input folder
    cd(inFolder);
    param;
    cd('../../../../../code/ETHZ/');
    disp(['get filters']);
    [allfilter, allsymbol] = makefilter(scale, norient);  % generate Gabor filters
    h = (size(allfilter{1}, 1)-1)/2;  % half size of Gabor
    C = corr(allfilter, epsilon);    % generate the inhibition maps
    
    I=[];Isize=[];
    disp(['get training images']);
   
    %% load negative data
    prepNegData;
    
    %% load positive data
    load([inFolder '/Imgs/posImg.mat']);
    I = IPos;
    Icolor = I;
    nimage  = numel(I);
    for i = 1:nimage
        I{i}=im2single(rgb2gray(Icolor{i}));
        Isize(i, :) = size(I{i});
    end
    sx = min(Isize(:, 1));
    sy = min(Isize(:, 2));
    Sx = sx*ones(1, nimage); 
    Sy = sy*ones(1, nimage);
    disp(['get S1 maps ']);
    for i = 1:nimage
        S1map = applyfilterfftsame(I(i), allfilter);  % filter training images
        ClocalNormalizeSingle(Sx(i),Sy(i),norient,h,localHalfx,localHalfy,...
            S1map,thresholdFactor);
        CsigmoidSingle(1,Sx,Sy,norient,saturation,S1map);
        outFileName = [workFolder  '/S1_' num2str(i,'%04d') 'ori' '.mat'];
        save(outFileName,'S1map');
    end
    
    disp(['get S1 maps and initialize M1 and S2 maps']);
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
    
    % Although the following code are done under an EM framwork, 
	% we did not use EM because we just run ONE interation and ONE category.
    %% begin the EM iteration
    scoreMat=zeros(nimage,numCluster);
    
    % random start;
    cluster = ceil(rand(nimage,1)*numCluster);
    clusterTemplate = cell(numCluster,1);

    for itrEM = 1:numItr
        Mstep; % learn the template, and adjust weight
        % Estep; % compute score on training images, it is not used currently.
    end
    selCurves = clusterTemplate{1};
    fileName = [outFolder '/' categories(indCate).name '_template.mat'];
    save(fileName,'selCurves','sx','sy','lambda','logZ',...
        'nAngle','nLength','Lrange2','Orange2','LRWeight');
    
end
