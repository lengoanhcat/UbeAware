% Main file for: Image Representation by Active Primitives
% To change the model parameters, edit the file "default_param.m"
% It is composed of both Mex-C routines and .m routines.

%% Compile Mex-C codes
mex ClocalNormalizeSingle.c % Local normalization
mex CMax1Single.c;  % compute the M1 map
mex CSum2Single.c;  % compute the S2 map
mex CBackM1andInhibition.c   % for inhibition and localization


inPath  = '../../data/singleImage';
outPath = '../../results/singleImage';

if ~exist(outPath) 
    mkdir(outPath);
end
% look up image files
imgList= dir([inPath '/0.jpg']);
numImage = length(imgList);

for iI= 1:numImage
    tic
    %% read image and parameters
    bike_param;
    disp(['loading image ' imgList(iI).name ])
    img_name = imgList(iI).name(1:end-4);
    
    img = imread([inPath '/' imgList(iI).name]);
    img = imresize(img,resizeFactor,'nearest');
    I= cell(1);
    I{1}=rgb2gray(img);
    [sx sy]=size(I{1});% image height and width

    
    %% build filters
    disp(['build filters']);
    [allfilter, allsymbol] = makefilter(scale, nOrient);  % generate Gabor filters
    h = (size(allfilter{1}, 1)-1)/2; % half size of Gabor
    C = corr(allfilter, epsilon);    % generate the inhibition maps
    
    nImage  = 1;
    %% Allocate memeory and compute S1 M1 S2
    S1map = cell(nOrient,1);
    M1map = cell(nOrient,1);
    S2map = cell(nLength,nAngle,nOrient);
    sym = zeros(sx,sy,'single'); % Image for symbolic prototypes
    Asym = cell(1);Asym{1}= zeros(sx,sy,'single'); % Image for deformed sketches
    
    
    for iO = 1 : nOrient
        M1map{iO} = NEGMAX+zeros(sx,sy,'single');% initialize M1maps
    end
    for iOri = 1:nOrient
        for iAng = 1: nAngle
            for iLen = 1:nLength
                S2map{iLen,iAng,iOri} = NEGMAX+ zeros(sx,sy,'single');% initialze S2maps
            end
        end
    end
    
    disp(['S1 maps...'])
    
    S1map= applyfilterfftsame(I(1), allfilter);  % filtering
    ClocalNormalizeSingle(sx,sy,nOrient,h,...
        localHalfx,localHalfy,S1map,thresholdFactor); % localNormalization
    sat = saturation;
    for iEl = 1:nOrient
        S1map{iEl}=sat.*(2./(1+exp(-2/sat.* S1map{iEl}))-1); % sigmoid transform
    end
    
    disp([num2str(toc) ' seconds up to now'])
    disp('M1...')
    CMax1Single(1,nOrient, h, S1map, M1map,Lrange, Orange, sx,sy);
    
    toc
    disp('S2...')
    CSum2Single(1, nOrient, nLength,nAngle, M1map, S2map, sx,sy, h,lambda,logZ);
    disp([num2str(toc) ' seconds up to now'])
    
    
    %% convert_index
    disp(['Convert index ...'])
    shiftS2Index;
    
    disp([num2str(toc) ' seconds up to now'])
    %% compute M2 maps
    disp(['M2 maps ...'])
    computeM2;
    disp([num2str(toc) ' seconds up to now'])
    
    %% S3PP map
    disp(['S3 maps...'])
    computeS3;
    disp([num2str(toc) 'seconds up to now'])
    
    %% non maximum surpression
    disp(['non-max supression...'])
    %X,Y: position of Corner
    %ORI_I, ORI_J: two arm orientations of a corner
    [X Y ORI_I ORI_J]=sparsify(S3PP,... %S3maps
        2*h+1,... %spatial neighorhood
        3,... %arm orientation neighborhood
        2.3*arcScoreTh,... % corner score threshold
        nOrient,sx,sy);% number of orienations, image size.
    disp([num2str(toc) ' seconds up to now'])
    nCorner = length(X); % number of corners
    
    
    
    %% trace back to S2 map
    % index of arcs used to compose corner
    arcX = zeros(2*nCorner,1);  % position
    arcY = zeros(2*nCorner,1);
    aOri = zeros(2*nCorner,1);  % orientation
    aLen = zeros(2*nCorner,1);  % length
    aAng = zeros(2*nCorner,1);  % curvature
    
    for i = 1:nCorner
        % find corresponding arcX,arcY,aOri,aLen,aAng
        % one arm
        iLen = ARG_M2PP_LEN{ORI_I(i)}(X(i),Y(i));
        iAng = ARG_M2PP_ANG{ORI_I(i)}(X(i),Y(i));
        map_vec  = backShiftMap{iLen,iAng,ORI_I(i)};
        aLen(i*2-1)=map_vec(1);
        aAng(i*2-1)=map_vec(2);
        aOri(i*2-1)=map_vec(3);
        dx =map_vec(4);
        dy =map_vec(5);
        arcX(i*2-1) = X(i)-dx;
        arcY(i*2-1) = Y(i)-dy;
        
        % the other arm
        iLen = ARG_M2PP_LEN{ORI_J(i)}(X(i),Y(i));
        iAng = ARG_M2PP_ANG{ORI_J(i)}(X(i),Y(i));
        map_vec  = backShiftMap{iLen,iAng,ORI_J(i)};
        aLen(i*2)=map_vec(1);
        aAng(i*2)=map_vec(2);
        aOri(i*2)=map_vec(3);
        dx =map_vec(4);
        dy =map_vec(5);
        arcX(i*2) = X(i)-dx;
        arcY(i*2) = Y(i)-dy;
    end% nCorner
    
    %% Trace back to S1map and draw the curves
    for iC = 1:nCorner
        for j =-1:0
            maxi =[aOri(iC*2+j) aAng(iC*2+j) aLen(iC*2+j) ...
                arcX(iC*2+j) arcY(iC*2+j)]; % vectorized index of arc
            maxi(1:3) = maxi(1:3)-1; % adapt to index system of C language
            imaxi = maxi; % deformed arc index, meaningless here.
            ds_respVec = zeros(2*imaxi(3)+1,1,'single'); % h(r) of each Gabor element
            
            CBackM1andInhibition(1, nOrient, nLength,nAngle,S1map, C,...
                sx,sy, h, Lrange, Orange,...
                maxi, imaxi, ds_respVec, sym, Asym, allsymbol(1, :),0);
            
        end%j
        % this time, inhibition
        for j =-1:0
            maxi =[aOri(iC*2+j) aAng(iC*2+j) aLen(iC*2+j) arcX(iC*2+j) arcY(iC*2+j)];
            maxi(1:3) = maxi(1:3)-1;
            imaxi = maxi;
            ds_respVec = zeros(2*imaxi(3)+1,1,'single');
            
            
            CBackM1andInhibition(1, nOrient, nLength,nAngle,S1map, C,...
                sx,sy, h, Lrange, Orange,...
                maxi, imaxi, ds_respVec, sym, Asym, allsymbol(1, :),1);
            
        end%j
    end% iC
    % copy corner lelvel symbolic representation
    cnSym = sym + 0.0;
    cnAsym = Asym{1}+ 0.0;
    
    
    %% begin arc pursuit
    disp('Selected Arcs:')
    disp('orient,       curvature,      length,         X,        Y,      score')
    selArcs = [];
    while true
        % M1,S2,
        CMax1Single(1,nOrient, h, S1map, M1map,Lrange, Orange, sx,sy);
        
        CSum2Single(1, nOrient, nLength,nAngle, M1map, S2map, sx,sy, h,lambda,logZ);
        
        % select best arc
        pooledMax = zeros(size(S2map)); % max scores on each single map
        pooledMaxInd = zeros(size(S2map)); % arg max of each single map
        for iEl = 1:numel(S2map)
            [val ind]=max(S2map{iEl}(:));
            pooledMax(iEl) = val;
            pooledMaxInd(iEl)=ind;
        end
        [val ind]=max(pooledMax(:));
        if val< arcScoreTh
            break;
        end
        
        % trace back [length,curvature,orient]
        [maxi(3),maxi(2),maxi(1)] = ind2sub(size(S2map),ind);
        maxi = maxi-1;
        ind = pooledMaxInd(maxi(3)+1,maxi(2)+1,maxi(1)+1);
        [maxi(4),maxi(5)]=ind2sub([sx,sy],ind);
        maxi(6) = val;
        selArcs = [selArcs; maxi];
        disp([ num2str(maxi)]);
        
        % inhibition
        imaxi = maxi;
        ds_respVec = zeros(2*imaxi(3)+1,1,'single');
        CBackM1andInhibition(1, nOrient, nLength,nAngle,S1map, C,...
            sx,sy, h, Lrange, Orange,...
            maxi, imaxi, ds_respVec, sym, Asym, allsymbol(1, :),1);
    end
    disp([ 'totally ' num2str(toc) ' seconds for current image']);
    
    disp('visualization...')
    % outpute original image
    imwrite(img,[outPath '/' imgList(iI).name(1:end-4)  '_0_org.png']);
    %imwrite(rgb2gray(img),[outPath '/' imgList(iI).name(1:end-4)  '_1_gray.png']);

    % corner saliency map
    scoreImg = max(0,M3PP);
    scoreImg = uint8(255*scoreImg/max(M3PP(:)));
    imwrite(255-scoreImg,[outPath '/' imgList(iI).name(1:end-4)  '_2_corner_score.png']);
    %save([outPath '/' imgList(iI).name(1:end-4) 'M3PP.mat'],'M3PP');
    
    % draw corners on original image
    imshow(img);
    hold on
    for i = 1:nCorner
        len = 15;
        alpha = (ORI_I(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'r.-','linewidth',3);
        
        alpha = (ORI_J(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'g.-','linewidth',3);
        
    end
    %set(gcf,'PaperPositionMode','auto');
    %saveas(gcf,[outPath '/' imgList(iI).name(1:end-4)  '_3_corner.png']);
    close
    
    % corners on blank image
    figure;
    imshow(ones(sx,sy));
    hold on
    line([1, 1,sy,sy,1],[1,sx,sx,1,1],'Color','k')
    for i = 1:nCorner
        len = 15;
        alpha = (ORI_I(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'r.-','linewidth',3);
        
        alpha = (ORI_J(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'g.-','linewidth',3);
        
    end
    set(gcf,'PaperPositionMode','auto');
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_8_corner_prototype.png'];
    saveas(gcf,filename);
    close;
    
    
    % prototype symbols
    imshow(1-cnSym,[0 1]);
    hold on
    for i = 1:nCorner
        len = 15;
        alpha = (ORI_I(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'r.-','linewidth',3);
        
        alpha = (ORI_J(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'g.-','linewidth',3);
        
    end
    set(gcf,'PaperPositionMode','auto');
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_4_corner_prototype.png'];
    saveas(gcf,filename);
    close;
    
    % deformed sketches
    imshow(1-cnAsym,[0 1]);
    hold on
    for i = 1:nCorner
        % visuallize corner points
        plot([Y(i)], [X(i)],'ro','linewidth',4);
    end
    set(gcf,'PaperPositionMode','auto');
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_6_corner_sym.png'];
    saveas(gcf,filename);
    
    
    
    % all active primitive prototypes
    imshow(1-sym,[0 1]);
    hold on
    for i = 1:nCorner
        len = 15;
        alpha = (ORI_I(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'r.-','linewidth',3);
        
        alpha = (ORI_J(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'g.-','linewidth',3);
        
    end
    set(gcf,'PaperPositionMode','auto');
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_5_all_prototype.png'];
    saveas(gcf,filename);
    close;
    
    
    % deformed version of all sketches
    imshow(1-Asym{1},[0 1]);
    hold on
    for i = 1:nCorner
        % visuallize corner points
        plot([Y(i)], [X(i)],'ro','linewidth',3);
    end
    set(gcf,'PaperPositionMode','auto');
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_7_all_sym.png'];
    saveas(gcf,filename);
    
    %filename = [outPath '/' imgList(iI).name(1:end-4)  '_7_all_sym.pdf'];
    %saveas(gcf,filename);
    
    
    filename = [outPath '/' imgList(iI).name(1:end-4)  '_statement.txt'];
    pFile = fopen(filename,'w');
    fprintf(pFile,'nCorner:%d,nArc:%d',nCorner,size(selArcs,1));
    fclose(pFile);
    close all;
end
