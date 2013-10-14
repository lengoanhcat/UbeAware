%function [S2PPmap backTrackMap]=shiftS2Index (S2map,nAngle,nOrient,nLength,h)
% create Gabor positions for each arc
startDx = zeros(nAngle,nOrient,nLength);
startDy = zeros(nAngle,nOrient,nLength);
endDx  = zeros(nAngle,nOrient,nLength);
endDy  = zeros(nAngle,nOrient,nLength);
startO =zeros(nAngle,nOrient,nLength);
endO   =zeros(nAngle,nOrient,nLength);

% link from endPoint index to center index.
backTrackMap = cell(nLength,nAngle,nOrient);

% following loop is do the same thing as arcShift function in Mex code
PI = 3.1415926;
cAngleLimit = (nAngle-1)/2;
for iOri = 0:nOrient-1  % mimic the C language index
    for iAng = 0:nAngle-1
        dAngle = -cAngleLimit+iAng;
        startDx(iAng+1,iOri+1,1)=0;
        startDy(iAng+1,iOri+1,1)=0;
        endDx(iAng+1,iOri+1,1)=0;
        endDy(iAng+1,iOri+1,1)=0;
        startO(iAng+1,iOri+1,1)=iOri+1;
        endO(iAng+1,iOri+1,1)=iOri+1;
        for l = 1:nLength-1
            startO(iAng+1,iOri+1,l+1)= (iOri+l*dAngle)+1;
            if(startO(iAng+1,iOri+1,l+1) < 1)
                startO(iAng+1,iOri+1,l+1)=startO(iAng+1,iOri+1,l+1)+nOrient;
            elseif startO(iAng+1,iOri+1,l+1) > nOrient;
                startO(iAng+1,iOri+1,l+1)=startO(iAng+1,iOri+1,l+1)-nOrient;
            end
            theta = PI/nOrient*(iOri+l*dAngle);  % no need to change to iOri+1
            dTheta = dAngle*PI/nOrient;
            startDx(iAng+1,iOri+1,l+1)=startDx(iAng+1,iOri+1,l)-round((h*1.8)*cos(dTheta/2)*sin(theta-dTheta/2));
            startDy(iAng+1,iOri+1,l+1)=startDy(iAng+1,iOri+1,l)+round((h*1.8)*cos(dTheta/2)*cos(theta-dTheta/2));
            
            endO(iAng+1,iOri+1,l+1) = (iOri-l*dAngle)+1;
            if endO(iAng+1,iOri+1,l+1) < 1
                endO(iAng+1,iOri+1,l+1)=endO(iAng+1,iOri+1,l+1)+nOrient;
            elseif endO(iAng+1,iOri+1,l+1) > nOrient;
                endO(iAng+1,iOri+1,l+1)=endO(iAng+1,iOri+1,l+1)-nOrient;
            end
            theta = PI/nOrient*(iOri+nOrient-l*dAngle);
            endDx(iAng+1,iOri+1,l+1)=endDx(iAng+1,iOri+1,l)-round((h*1.8)*cos(dTheta/2)*sin(theta+dTheta/2));
            endDy(iAng+1,iOri+1,l+1)=endDy(iAng+1,iOri+1,l)+round((h*1.8)*cos(dTheta/2)*cos(theta+dTheta/2));
        end
    end
end

%% begin convering the index
S2PPmap = cell(nLength,nAngle,2*nOrient); % S2map, using end point index
midAngle = (nAngle+1)/2;
for iAngle= midAngle-dCurvature:midAngle+dCurvature
    dAngle =  iAngle-(nAngle+1)/2;
    for iOri = 1:nOrient        
        for iLength = minLength:nLength %only convert arcs longer than mininum length
            % one end point
            dx = startDx(iAngle,iOri,iLength);
            dy = startDy(iAngle,iOri,iLength);
            newO = iOri + (iLength-1)*dAngle+nOrient;
            if newO <1
                newO = newO + 2* nOrient;
            elseif newO >2*nOrient
                newO = newO -2*nOrient;
            end
            newAngle= nAngle-(iAngle-1);
            
            x_min = max(1,1+dx);y_min=max(1,1+dy);
            x_max = min(sx,sx+dx); y_max= min(sy,sy+dy);
            S2PPmap{iLength,newAngle,newO} = -1e10+zeros(sx,sy,'single');
            S2PPmap{iLength,newAngle,newO}(x_min:x_max,y_min:y_max)= ...
                S2map{iLength,iAngle,iOri}(x_min-dx:x_max-dx,y_min-dy:y_max-dy);
            backTrackMap{iLength,newAngle,newO}=[iLength,iAngle,iOri,dx,dy];
            
            
            % the other end point
            dx = endDx(iAngle,iOri,iLength);
            dy = endDy(iAngle,iOri,iLength);
            newO = iOri-(iLength-1)*dAngle;
            if newO <1
                newO = newO + 2* nOrient;
            elseif newO >2*nOrient
                newO = newO -2*nOrient;
            end
            x_min = max(1,1+dx); y_min=max(1,1+dy);
            x_max = min(sx,sx+dx);y_max= min(sy,sy+dy);
            S2PPmap{iLength,iAngle,newO} = -1e10 + zeros(sx,sy,'single');
            S2PPmap{iLength,iAngle,newO}(x_min:x_max,y_min:y_max)= ...
                S2map{iLength,iAngle,iOri}(x_min-dx:x_max-dx,y_min-dy:y_max-dy);
            backTrackMap{iLength,iAngle,newO}=[iLength,iAngle,iOri,dx,dy];
            
        end
    end
end
clear startDx startDy startO endDx endDy endO;
backShiftMap = backTrackMap;
