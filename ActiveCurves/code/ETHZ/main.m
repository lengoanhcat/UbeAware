for run_id = 1:1
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',run_id+4));
    
    out_path = ['../../data/cache/ETHZ/extractData' num2str(run_id)];
    buffer_path = ['../buffer' num2str(run_id)];
    global FaceDB; FaceDB = 0;

    % arguments in cropData are:
    % Source image directory
    % Image output directory
    % Window height
    % Bounadary padding size
    % Proportion of images for training
    % Number of positive images for each negative image
%     disp('cropping Bottles ...')    
%     cropDataRandom('../../data/ETHZ/Bottles',[ out_path '/Bottles'],130,10,0.5,0.5);
%     % copy the parameter file
%     copyfile('Bottles_param.m',[ out_path '/Bottles/param.m'])    
%     
%     disp('cropping Mugs ...')
%     cropDataRandom('../../data/ETHZ/Mugs',[ out_path '/Mugs'],85,10,0.5,0.5);
%     copyfile('./Mugs_param.m',[ out_path '/Mugs/param.m']);
%         
%     disp('cropping Giraffes ...')
%     cropDataRandom('../../data/ETHZ/Giraffes',[ out_path '/Giraffes'],145,10,0.5,0.5);
%     copyfile('./Giraffes_param.m',[ out_path '/Giraffes/param.m']);
%     
%     disp('cropping Applelogos ...')
%     cropDataRandom('../../data/ETHZ/Applelogos',[ out_path '/Applelogos'],115,10,0.5,0.5);
%     copyfile('./Applelogos_param.m',[ out_path '/Applelogos/param.m']);
%     
%     disp('cropping Swans ...')
%     cropData('../../data/ETHZ/Swans',[ out_path '/Swans'],55,10,0.5,0.5)
%     copyfile('./Swans_param.m',[ out_path '/Swans/param.m']);
%     
%     FaceDB = 1;
%     disp('cropping Faces ...')
%     cropData('../../data/ETHZ/Faces',[ out_path '/Faces'],55,10,0.5,0.5);
%     copyfile('./Faces_param.m',[ out_path '/Faces/param.m']); 
    
    FaceDB = 1;
    disp('cropping Faces ...')
%     cropData('../../data/ETHZ/FrontalFaces3',[ out_path '/FrontalFaces3'],55,0,0.2,0.8);
%     loadData; % temporarily load train and test data
    copyfile('./emotion_param.m',[ out_path '/disgust/param.m']);
    % learn template and adjust feature weight
    ABline;
    
    % performe detection task.
    % You can run each of the following five lines in a separate matlab thread.
%     testCategoryV2('Bottles',buffer_path,run_id);
%     testCategoryV2('Mugs',buffer_path,run_id);
%     testCategoryV2('Giraffes',buffer_path,run_id);
%     testCategoryV2('Applelogos',buffer_path,run_id);
%     testCategoryV2('Swans',buffer_path,run_id);
%     testCategoryV2('Faces',buffer_path,run_id);
%     testCategoryV2('disgust',buffer_path,run_id);
end
% produce the curves presented in paper
% finalCurve;