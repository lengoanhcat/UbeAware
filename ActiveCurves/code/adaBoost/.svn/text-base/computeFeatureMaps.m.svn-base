% This m file compute the features on M1 maps, M2^+ maps, and M3 maps.
% After computing the feature maps, these maps are concantated into a vector,
% There boundaries for each map are saved in variable bound.
% These features will be selected by adaBoost to form models for object classification. 

if ~exist(cachePath)
    mkdir(cachePath);
end
if ~exist(outPath)
    mkdir(outPath);
end


%% read image and parameters
default_param;
label = {'pos','neg'};
%% build filters
disp(['build filters']);
[allfilter, allsymbol] = makefilter(scale, nOrient);  % generate Gabor filters
h = (size(allfilter{1}, 1)-1)/2; % half size of Gabor
C = corr(allfilter, epsilon);    % generate the inhibition maps


for indCate = 1:length(cl);
    imgList= dir([inPath '/' cl{indCate} '/*.jpg']);
    numImage = length(imgList);
    for iI= 1:numImage
	dataFileName = [cachePath '/' label{indCate} '_' num2str(iI,'%04d') '.mat'];
	if exist(dataFileName,'file')
		continue;
	end
	save(dataFileName,'iI');
        tic
        disp(['loading image ' imgList(iI).name ])
        img_name = imgList(iI).name(1:end-4);
        
        img = imread([inPath '/' cl{indCate} '/' imgList(iI).name]);
        img = imresize(img,[height width],'nearest');
        I= cell(1);
        if size(img,3)==3
            img = rgb2gray(img);
        end
        I{1}=img;
        [sx sy]=size(I{1});% image height and width
        
        
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
        sat = single(saturation);
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
        
		toc
		disp('M2...');
		M2map = cell(nLength,nAngle,nOrient);
		for imap = 1:nLength*nAngle*nOrient
			M2map{imap}=-1e10+zeros(sx,sy,'single');			
		end
        CMax2Single(1, nOrient, h, nLength, nAngle, ...
        		S2map, M2map, Lrange2, Orange2, Crange2,LenRange2,...
            	sx,sy);


        %% convert_index
        disp(['Convert index ...'])
	    shiftS2Index;
        
        disp([num2str(toc) ' seconds up to now'])
        %% compute M2 maps
        disp(['M2PP maps ...'])
        computeM2;
        disp([num2str(toc) ' seconds up to now'])
        
        %% S3PP map
        disp(['S3 maps...'])
        computeS3;
        disp([num2str(toc) 'seconds up to now'])
        
        %% lump the maps into a sparse vector;
        bound = [];
        %m1Vec = sparse(double(reshape(cell2mat(M1map),1, [])));
        m1Vec = reshape(cell2mat(M1map),1, []);
        bound = [bound; length(m1Vec)];
        
        %m2Vec = sparse(double(reshape(cell2mat(M2PPmap),1,[])));
        m2Vec = reshape(cell2mat(M2map),1,[]);
        bound = [bound; bound(end)+length(m2Vec)];
        
        m3Vec = [];
	    m3VecBound = [];
	    n_hood = ones(2*dxCorner+1,2*dyCorner+1);
        for iOri = 1:2*nOrient % orientation for first arm
            for j=iOri+minAng:iOri+maxAng % possible orientation for second arm
                if j<1  % mod the orienation index to range [0 2K-1];
                    jOri = j+2*nOrient;
                elseif j>2*nOrient
                    jOri = j-2*nOrient;
                else
                    jOri = j;
                end
                %s3Vec = [s3Vec sparse(double(reshape(S3PP{iOri,jOri},1,[])))]; % compute S3Map
		m3VecBound = [m3VecBound; [iOri,jOri,length(m3Vec)]];
		m3map_i_j = imdilate(S3PP{iOri,jOri},n_hood);
                m3Vec = [m3Vec reshape(m3map_i_j,1,[])]; % compute S3Map
            end
        end
        bound = [bound; bound(end)+length(m3Vec)];
        
        fVec =[m1Vec m2Vec m3Vec];
        save(dataFileName, 'I','img_name','sx','sy','fVec','bound','m3VecBound');
        toc
	    fprintf('\n\n')
        
    end% for image
end % for category
