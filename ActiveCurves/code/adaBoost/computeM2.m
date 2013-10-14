%% Routines to compute M2 maps
M2PPmap = cell(2*nOrient,1); % Max2 maps
ARG_M2PP_ANG=cell(2*nOrient,1); % Arg-max maps, for curvature
ARG_M2PP_LEN=cell(2*nOrient,1); % Arg-max maps, for arc length
for iEl = 1:2*nOrient % ininalize maps
    M2PPmap{iEl}=-1e10+zeros(sx,sy,'single');
    ARG_M2PP_ANG{iEl}=zeros(sx,sy,'uint8');
    ARG_M2PP_LEN{iEl}=zeros(sx,sy,'uint8');
end

for iOri = 1:2*nOrient
    m_map = M2PPmap{iOri}; % current M2 map
    arg_ang = ARG_M2PP_ANG{iOri}; % current arg-max map, for curvature
    arg_len = ARG_M2PP_LEN{iOri}; % current arg-max map, for arc length
    midAngle = (nAngle+1)/2; % index corresponding to \rho = 0;
    for iAngle = midAngle-dCurvature:midAngle+dCurvature
        for iLength = minLength:nLength
            m_ind = (m_map<S2PPmap{iLength,iAngle,iOri});
            m_map(m_ind) = S2PPmap{iLength,iAngle,iOri}(m_ind);
            arg_ang(m_ind)=iAngle;
            arg_len(m_ind)=iLength;
        end
    end
    M2PPmap{iOri}=m_map;
    ARG_M2PP_ANG{iOri}=arg_ang;
    ARG_M2PP_LEN{iOri}=arg_len;
end
clear S2PPmap; % clear end point indexed S2 maps
