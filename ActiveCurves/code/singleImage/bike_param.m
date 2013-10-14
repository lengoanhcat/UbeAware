scale = .8;       % s; scale of Gabors, length = 17 pixels                        
nOrient = 25;     % K; number of orientations              
Lrange = 1;       % b_4; location shift range of Gabor Primitives                  
Orange = 1;       % b_5; orientation deform range of Gabor Primitives
nAngle = 9;       % \rho; number of different curvatures(should be odd number)  
nLength = 4;      % l; number of length of line segments  
dCurvature = 1;   % b_1; because matlab index start from 1
minLength = 2;    % b_2+1; because the range is inclusive, and indexed from 1             
minAng = 8;       % b_3; mininum corner angle, inclusive                                                            
maxAng = 15;      % b_3; maximum corner angle, inclusive 
lambda = 2.5;     % \lambda; expected model coefficient                        
logZ = 10.4;      % \log Z; this is computed offline                                
arcScoreTh=3*(5.6423*lambda-logZ); % T; threshold on arc score
cnScoreTh = 2.3*arcScoreTh;        % b_6; corner score threshold

%% Parameters for local normalization
localHalfx=round(40*scale); % halflength of local normalization area
localHalfy=round(40*scale); % halflength of local normalization area
Upperbound = 16.; % saturation value
thresholdFactor = 0.01;
%% Misc parameters
saturation = 6.0; % saturation value for sigmoid transformation
epsilon = .01;    % allowed correlation between selected Gabors
NEGMAX=-1e10;     % a value for minnimum value
resizeFactor = 0.8; % used to resize the image, if necessary