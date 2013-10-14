S3PP = cell(2*nOrient,2*nOrient); %S3 maps
M3PP = zeros(sx,sy); % M3 maps, for visualization
for iOri = 1:2*nOrient % orientation for first arm
    for j=iOri+minAng:iOri+maxAng % possible orientation for second arm
        if j<1  % mod the orienation index to range [0 2K-1];
            jOri = j+2*nOrient;
        elseif j>2*nOrient
            jOri = j-2*nOrient;
        else
            jOri = j;
        end
        S3PP{iOri,jOri}=M2PPmap{iOri}+M2PPmap{jOri}; % compute S3Map
        m_ind = S3PP{iOri,jOri}>M3PP;
        M3PP(m_ind)=S3PP{iOri,jOri}(m_ind);
    end
end
