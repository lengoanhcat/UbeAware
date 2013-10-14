for iCluster = 1:numCluster
    indImage = find(cluster==iCluster); % find images belong to this cluster
    posRespMap =[];negRespMap=[];

    disp(['Mstep, iter ' num2str(itrEM), ' cluster:' num2str(iCluster) ])
    sym = zeros(sx, sy,'single');
    Asym = cell(1, nimage);   % symbolic plot for each image with active Gabor
    for (i = 1 : nimage)
        Asym{i} = zeros(sx, sy,'single');
    end

    for i = 1 : length(indImage)
        % refresh the S1map
        sourceFile = [workFolder  '/S1_'  num2str(indImage(i) ,'%04d') 'ori' '.mat'];
        destFile = [workFolder  '/S1_'  num2str(indImage(i),'%04d')  '.mat'];
        copyfile(sourceFile,destFile);
    end


    sk=0; % number of selected sketches
    maxi = zeros(1, 6); imaxi = zeros(1, 6);
    while(sk <=Totalsketch)
        %% refresh the pooled M2
        pooledMap = cell(nLength,nAngle,norient); % shared M2^+ maps
        pooledMax = zeros(nLength,nAngle,norient); % maxium response on each M2^+ map
        pooledMaxInd = zeros(nLength,nAngle,norient); % location of maximum response on each M2^+ map
        for iEl = 1:nLength*nAngle*norient
            pooledMap{iEl} = zeros(sx,sy,'single'); % initialize shared M2^+ map
        end

		% compute M2^+ maps for each image, and add them to shared M2^+ maps
        for i = 1 : length(indImage)
            disp([ 'to ' num2str(indImage(i)) 'th image'])
            tic
            inFileName = [workFolder  '/S1_'  num2str(indImage(i),'%04d')  '.mat'];
            load(inFileName);

            CMax1Single(1,norient, S1map, M1map, Lrange, Orange, Sx, Sy);

            CSum2Single(1, norient, nLength,(nAngle-1)/2, M1map, S2map, sx, sy, h,lambda,logZ);

            CMax2Single(1, norient, nLength, nAngle, S2map, M2map, Lrange2, Orange2, Sx, Sy);
            %%update the pooled M2map
            for iEl = 1:nLength*nAngle*norient
                pooledMap{iEl}=pooledMap{iEl}+M2map{iEl};
            end
            disp( ['totally ' num2str(toc) ' seconds for current image']);

        end
        disp(['get the next sketch']);
        curveSelection;

        disp(['retrive arg-max and inhibit S1 map']);
        respVects = zeros(2*maxi(3)+1,length(indImage),'single');
        respVect = zeros(2*maxi(3)+1,1,'single');
        for i = 1 : length(indImage)
            disp([ 'to ' num2str(indImage(i)) 'th image'])

            inFileName = [workFolder  '/S1_'  num2str(indImage(i),'%04d')  '.mat'];
            load(inFileName);
            CMax1Single(1,norient, S1map, M1map, Lrange, Orange, Sx, Sy);
            CSum2Single(1, norient, nLength,(nAngle-1)/2, M1map, S2map, sx, sy, h,lambda,logZ);

			
            CBackMax2Single(1, norient, nLength, nAngle, S2map, h, Lrange2, Orange2, sub, Sx, Sy, maxi, imaxi);
            CBackM1andInhibition(1, norient, nLength, nAngle, S1map, C, sx, sy, h, Lrange, Orange,...
                maxi, imaxi, respVect, sym, Asym(indImage(i)), allsymbol(1, :));
			% save inhibited S1 maps into disk.			            
			outFileName = [workFolder  '/S1_'  num2str(indImage(i),'%04d')  '.mat'];
            save(outFileName,'S1map');
            respVects(:,i) = respVect;
        end
        posRespMap = [posRespMap respVects'];

        sk = sk + 2*maxi(3)+1;
        selCurves=[selCurves;maxi];

    end
    clusterTemplate{iCluster}=selCurves;

    disp(['retrive responses on negative data']);
    negRespMap=[];
	for iArc = 1:size(selCurves,1)
        maxi=selCurves(iArc,:);
        respVect = zeros(2*maxi(3)+1,1,'single');
        respVects = zeros(2*maxi(3)+1,nNeg,'single');
        negSym = zeros(size(sym),'single');
        negAsym = cell(1);
        negAsym{1}=zeros(size(sym),'single');
        for iNeg = 1:nNeg
            outFileName = [workFolder  '/S1_neg_' num2str(iNeg,'%04d') 'ori' '.mat'];
            load(outFileName);
            
            outFileName = [workFolder  '/S2_neg_' num2str(iNeg,'%04d') 'ori' '.mat'];
            load(outFileName);

            CBackMax2Single(1, norient, nLength, nAngle, S2map, h, Lrange2, Orange2, sub, Sx, Sy, maxi, imaxi);
            CBackM1andInhibition(1, norient, nLength, nAngle, S1map, C, sx, sy, h, Lrange, Orange,...
                maxi, imaxi, respVect, negSym, negAsym(1), allsymbol(1, :));
            respVects(:,iNeg) = respVect;
        end
        negRespMap = [negRespMap respVects'];
   end



    %% adjust weight using logistic regression
    param_s = 0;
    param_B = 1;
    param_c = 0.01;
    weight_of_pos = 1; %double(sum(y<0))/double(sum(y>0));
    % -s 0 is for L2-logistic regression, -B 1 is for bias term, -c is for tuning parameterples
    para=sprintf('-s %d -B %d -c %f -w1 %f -q',param_s,param_B,param_c,weight_of_pos); 

    Xt = [posRespMap;negRespMap];
    X = sparse(double(Xt));
    Y = ones( size(X,1), 1 );
    Y(1+size(posRespMap,1):end) = -1;    % negative examples
    model = lrtrain(Y,X,para);
    LRWeight=model.w;

    % OUTPUT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = length(indImage); bx =0; by = 0; cc = 3;
    col = 10; row = floor(K/col)+1;
    Iout = zeros(row*sx+(row-1)*(2*bx+cc), 2*col*sy+(col-1)*(2*by+cc))+255;
    towrite =  -sym;
    towrite = (255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
    Iout(1:sx, 1:sy) = towrite;
    i = -1;
    for (r = 1 : row)
        for (c = 1 : col)
            i = i+1;
            if ((i>=1)&&(i<=K))
                towrite = I{indImage(i)};
                towrite1 = (255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
                Iout((r-1)*sx+1+(r-1)*(2*bx+cc): r*sx+(r-1)*(2*bx+cc), (c-1)*2*sy+1+(c-1)*(2*by+cc): (c-1)*2*sy+sy+(c-1)*(2*by+cc)) = towrite1;
                towrite =  -Asym{indImage(i)};
                towrite2 = (255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
                %towrite2 = 0.5*towrite1 + 0.5*towrite2;
                Iout((r-1)*sx+1+(r-1)*(2*bx+cc) : r*sx+(r-1)*(2*bx+cc), (c-1)*2*sy+sy+(c-1)*(2*by+cc)+1 : (c-1)*2*sy+sy+(c-1)*(2*by+cc)+sy) =towrite2;
            end
        end
    end
    for (r = 1 : row-1)
        Iout(r*sx+(r-1)*(2*bx+cc)+bx+1:r*sx+(r-1)*(2*bx+cc)+bx+cc, :) = 0;
    end
    for (c = 1 : col-1)
        Iout(:, c*(sy*2)+(c-1)*(2*by+cc)+by+1:c*(sy*2)+(c-1)*(2*by+cc)+by+cc) = 0;
    end

    Iout = uint8(255 * (Iout-min(Iout(:)))/(max(Iout(:))-min(Iout(:))));

    fileName = [outFolder '/' categories(indCate).name 'itr' num2str(itrEM) 'cl' num2str(iCluster) '.png']
    imwrite(Iout,fileName);
    fileName = [outFolder '/' categories(indCate).name 'itr' num2str(itrEM) 'cl' num2str(iCluster) '_template.mat']
    save(fileName,'selCurves','sx','sy','lambda','logZ');
    fileName = [outFolder '/' categories(indCate).name 'itr' num2str(itrEM) 'cl' num2str(iCluster) '_all.mat']
    save(fileName);
    close;
end
