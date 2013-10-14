function [X Y ORI_I ORI_J]=sparsify(S3PP,r,ang,th,nOrient,sx,sy)
%Perform the operations described in footnote 1
%S3PP: S3 maps
%r: spatial neighborhood range
%ang: orientation neighborhood range
%th: score threshold

accX=[]; accY=[]; accI=[]; accJ=[];accVal = [];
thMask = ones(sx,sy,'single')*th;

for iOri = 1:2*nOrient
    for jOri = 1:2*nOrient
        if isempty(S3PP{iOri,jOri})
            continue;
        end
        map = S3PP{iOri,jOri};
        [I,J]=find(S3PP{iOri,jOri}>th);
        if isempty(I)
            continue;
        end
        K = sub2ind([sx,sy],I,J);
        map =sparse(I,J,double(map(K)));% sparse maps, only points pass the threshold 
        [x y]= nonmaxsuppts(map,1,th); % non-max supression, usnig peters code, see citation [11]
        ind = sub2ind([sx sy],x,y);
        val = S3PP{iOri,jOri}(ind);
        ori_j = ones(size(x))*jOri;
        ori_i = ones(size(x))*iOri;
        accX = [accX; x];
        accY = [accY; y];
        accI = [accI; ori_i];
        accJ = [accJ; ori_j];
        accVal = [accVal; val];
    end
end


% Another level of non-max supression
candMat = [accX accY accI accJ accVal];
selMat = [];
X=[];Y=[];ORI_I=[];ORI_J=[];
nPt = length(accX);
for iPt = 1:nPt
    % spatial neighrhood
    x = accX(iPt); y = accY(iPt);
    dist = sqrt((accX-x).^2+(accY-y).^2);
    spatial_neighbor = (dist<r);
    
    % orientation neighorhood
    ori_i = accI(iPt); ori_j = accJ(iPt);
    distI = min(abs(accI-ori_i),abs(accI+2*nOrient-ori_i));
    distI = min(distI,abs(ori_i+2*nOrient-accI));
    
    distJ = min(abs(accJ-ori_j),abs(accJ+2*nOrient-ori_j));
    distJ = min(distJ,abs(ori_j+2*nOrient-accJ));
    
    orient_neighbor = (distI<=ang)&(distJ<=ang);
    
    neighbor = orient_neighbor&spatial_neighbor;
    neighbor(iPt)=false;% exclude itself
    
    % see if it is maximum
    if sum(neighbor)>0
        val = accVal(iPt);
        if val<max(accVal(neighbor))
            continue;
        end
        % break ties 
        %if length(find(accVal(neighbor)==val))>0
        %    accVal(iPt)=-1e10;
        %    continue;
        %end
    end
    %push to result
    selMat = [selMat; candMat(iPt,:)];
end
% output
if ~isempty(selMat)
    X = selMat(:,1);Y=selMat(:,2);ORI_I=selMat(:,3);ORI_J=selMat(:,4);
end
