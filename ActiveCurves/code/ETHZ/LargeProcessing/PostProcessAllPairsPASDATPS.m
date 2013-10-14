function [S, Hvals] = PostProcessAllPairsPASDATPS(model_dir, test_dir, process_dir, postproc_pars, draw_pars, subset, score_fct, fullres_params)

% Reloads process files for all pairs of model and test images
% and carries out some post-processing
% (e.g export figures, count correct detections and false-pos)
%
% Returns statistics S. Useful for building them incrementally from several test directories.
%
% Input:
% - subset: which test images to post-process;
%           if 'non-trn'          -> process all images in test_dir but those in model.trn_ixs
%           if not given or 'all' -> process all images in test_dir
% - if fullres_params given ->
%   recompute hyp cues at full resolution (for this needs the params of shape matcher that produced them)
%
% Output:
% - Hvals = list of all hypotheses from all test images with the values of their cues, etc.
% - if score_fct given -> use it to re-score hypotheses; if not given -> keep scores as in the existing process files
%
% If model_dir is a cell-array containing 2 strings ->
% model_dir{1} = the model directory
% model_dir{2} = the class name (used to reconstruct the groundtruth filename),
%                using this option is useful when the model filenames doesn't reveal trivially the class name
%
%


% parse arguments
if nargin < 8
  fullres_params = false;
end
if nargin < 7
  score_fct = false;
end
if nargin < 5
  draw_pars = false;
end

% save current working directory
org_dir = pwd;

% resolve model_dir
if iscell(model_dir)
  class_name = model_dir{2};
  model_dir = model_dir{1};
else
  class_name = false;                       % will resolve the class name from the model filenames
end

% get list of model images
cd(model_dir);
models = dir('*.mat');
disp('Models found: ');
dir('*.mat');                               % must display this way, because variable is a struct
newline;

% get list of test images
cd(org_dir); cd(test_dir);
tests = dir('*.mat');
disp('Test images found: ');
dir('*.mat');             
newline;

subset_ixs = subset;                        % when user gives exact image indeces
if nargin < 5 || strcmp(subset,'all')
  subset_ixs = 1:length(tests);             % by default all images in test dir are processed
  subset = 'all';
end


% process all pairs
mix = 0;                                    % index of current model in 'models'
for m = models'
  mix = mix+1;
  m_base = m.name(1:(length(m.name)-4));
  S(mix) = InitStats;
  skipped = false;
  newline; newline;
  disp(['Processing model ' m_base]);
  disp('===========================================');
  cd(org_dir);
  model = load([model_dir '/' m.name]);
  if strcmp(subset,'non-trn')
    if isfield(model.obj, 'trn_ixs')
      subset_ixs = setdiff(1:length(tests), model.obj.trn_ixs);  % only process non-trn images
      disp('requested to process only non-training images');
      disp(['test subset = ' num2str(subset_ixs)]);
      newline;
    else
      error([mfilename ': requested to process only non-training images, but no .trn_ixs field found in model.']);
    end
  end

  Hvals = [];
  missed = zeros(1,length(subset_ixs));
  for test_img_ix = subset_ixs
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Loading
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    t = tests(test_img_ix);
    t_base = t.name(1:(length(t.name)-4));
    newline;
    disp(['Post-processing ' m_base ' in ' t_base]);
    disp('------------------------------------------------------------');
  
    % load data
    cd(org_dir);
    if not(islogical(draw_pars))
     disp('Loading test datafile');
     test = load([test_dir '/' t.name]); % if need to draw, need extra fields (.ifname, etc.) from test
    else
      % when not drawing -> just need test.obj.name -> avoid loading it -> BIG speedup
      disp('Not requested to draw anything -> speedup by not loading test mat-file.');
      test.obj.name = t_base;     
    end
   
    % load process-file
    process_fname = ['process_' m_base '_' t_base '.mat'];
    %process_fname = ['process_myentire1_' t_base '.mat'];  % ONLY FOR OLD JURIE-HORSES PROCESS FILES !
    cd(process_dir);
    if exist(process_fname, 'file')
      disp(['Loading ' process_fname]);
      clear trash;
      load(process_fname);
      if exist('trash','var')                               % if was a guard file -> skip post-processing
         disp(['Process file ' process_fname ' is a guard file --> skipping it']);
         continue;
      else
        disp(['successfully load process file ' process_fname]);
      end
    else
      disp(['Process file ' process_fname ' doesn''t exist --> skipping model ' m_base]);
      skipped = true;
      break;
    end
    
   
    %%%%%%%%%%%%%%%%%
    % Post-processing
    %%%%%%%%%%%%%%%%%
    if not(isempty(score_fct)) && not(islogical(score_fct))
      if not(islogical(fullres_params)) && not(isempty(fullres_params))
        disp('Computing full resolution cues');
        H = ComputeValsFullRes(1:length(H), H, model.obj, test.obj, fullres_params, true, true);
      end
      %
      % Hough-only mode
      if postproc_pars.hough_only
        disp('Hough-only mode');
        H = ConvertToHoughOnly(H, model.obj);
        clear temp; temp.H = H;
      else
        disp('Recomputing hypothesis scores');
        clear temp; temp.H = H;
        temp = RecomputeNrgsAndScores(temp, score_fct);
      end
      %
      % does all post-processing, including exporting figures and saving process file (if last param == true)
      cd(org_dir);
      temp = ConvertAllPASDATPS(model.obj, test.obj, temp, test_dir, process_dir, draw_pars, postproc_pars.remove_pars, true);
      H = temp.H;  C = temp.C;  dets = temp.dets;
    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute performance statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % partition detections in correct detections and false-positives
    % and compute scores of all positive images and negative images
    % for ROC plot for binary classification task
    cd(org_dir); cd(test_dir);
    if class_name
        gt_fname = [t_base '_' class_name '.groundtruth'];
    else  % resolve groundtruth filename from model filename m_base
        gt_fname = [t_base '_' m_base(3:(end-1)) 's.groundtruth'];    % trailing 's' for plural
    end
    %
    clear temp;
    [S(mix) temp corr_dets] = UpdateStats(S(mix), gt_fname, C, dets, 'all', postproc_pars.pascal);
    missed(test_img_ix) = temp;
    newline;

    % associate to each true-positive detection (decided by a BB intersection criterion)
    % a shape-matching accuracy (in [0,1]).
    if ~isempty(corr_dets) && postproc_pars.comp_accuracies                    % only do it if there's something to compute -> speedup
      cd(org_dir);
      [GTO gto_type] = LoadGTOutlines(test_dir, t_base, class_name);
      if not(isempty(GTO))
        corr_hixs = [];
        for corr_det = corr_dets
          corr_hixs = [corr_hixs C(corr_det).H(1)];
        end
        %ds = HypsDissim(H(corr_hixs), GTO, gto_type, model.obj, 'vx','GTO');  % meas = mean sym-chamf
        [ds ds_hyp ds_gto gto_ixs] = HypsDissim(H(corr_hixs), GTO, gto_type, model.obj, 'vx', 'GTO', 'perc_corr', 0.04);  % meas = perc corr loc pts
        S(mix).corr_acc = [S(mix).corr_acc ds];
        %S(mix).corr_acc = [S(mix).corr_acc ds_gto];
        disp(['accuracies:         ' num2str(ds)]);
        %disp(['accuracies:         ' num2str(ds_gto)]);
        if length(ds) < length(corr_dets)
          disp('Warning: less accuracy values than correct detections !');
        end
        % store which image and annotated obj within the image produced each corr_acc val
        S(mix).corr_ids = [S(mix).corr_ids [ones(1,length(gto_ixs))*test_img_ix; gto_ixs]];
        %
        % compute accuracy of ground-truth BBs
        Hgto_bb = MakeGTOBBShapes(GTO, gto_type, gto_ixs);
        [ds ds_hyp ds_gto gto_ixs] = HypsDissim(Hgto_bb, GTO, gto_type, false, 'vx', 'GTO', 'perc_corr', 0.04);
        disp(['accuracies of gt BB:         ' num2str(ds)]);
        S(mix).gt_acc = [S(mix).gt_acc ds];
      else
        disp('ground-truth outlines not found.');
        disp('skipping computation of accuracy of matched shapes.');
        keyboard;
        newline;
      end
    end

  end  % loop over test images



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % performance statistics for this model
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if not(skipped)
    disp(['Model ' m_base ' completed.']);
    if any(missed)
      disp('Missed instances: ');
      for tix = 1:length(subset_ixs)
        if missed(tix)
          disp([num2str(missed(tix)) ' missed in ' tests(subset_ixs(tix)).name(1:(length(tests(subset_ixs(tix)).name)-4))]);
        end
      end
    else
      disp('No instance missed !');
    end
  end

end % loop over models

% restore current working directory
cd(org_dir);
