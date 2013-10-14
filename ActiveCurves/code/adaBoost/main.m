% Main file for: Image Representation by Active Primitives
% To change the model parameters, edit the file "default_param.m"
% It is composed of both Mex-C routines and .m routines.

%% Compile Mex-C codes
mex ClocalNormalizeSingle.c % Local normalization
mex CMax1Single.c;  % compute the M1 map
mex CSum2Single.c;  % compute the S2 map
mex CBackM1andInhibition.c  % for inhibition and localization
mex CMax2Single.c;


%% Experiment for adaBoost on emotion disgust data
inPath  = '../../data/adaBoost';
outPath = '../../results/adaBoost/happiness';
cachePath = '../../data/cache/adaBoost/happiness';
cl = {'happiness','negativeImage'};
default_param; % load paramters
height = 172; width = 140; % define image size
computeFeatureMaps; % coompute the features and save them into disk
adaBoost_AB; % run adaBoost on active basis features
adaBoost_ARC; % run adaBoost on acitve arcs features;
adaBoost_ARC_CN; % run adaBoost on active arcs and active corners;
drawAUCComparison; % draw the comparison figure as in paper


% %% Experiment for adaBoost on leaves data
% inPath = '../../data/adaBoost';
% outPath = '../../results/adaBoost/leaves';
% cachePath = '../../data/cache/adaBoost/leaves';
% cl = {'horse','negativeImage'};% Folder name for positive and negative images
% default_param; % Load parameters
% height = 120; width = 180; % Define image size;
% computeFeatureMaps; % Compute the features and save them into disk
% adaBoost_AB; % Run adaBoost on active basis features
% adaBoost_ARC; % Run adaBoost on acitve arcs features;
% adaBoost_ARC_CN; % Run adaBoost on active arcs and active corners;
% drawAUCComparison; % draw the comparison figure as in paper