function DetectAllPairsPASDATPS(model_dir, test_dir, process_dir, det_pars, draw_pars, remove_pars, HGA, score_fct, skip_done, subset)

% Matches models in model_dir to test images in test_dir
% and saves output to process_dir.
%
% Each hypothesis energy is computed using 
% score_fct(1), and then converted to a score using score_fct(2)
%
% skip_done: do not process already done test images,
%            or those being done by another process at the same time
%
% subset: which test images to process;
%         if 'non-trn' -> all test images, but those in model.trn_ixs
%         if not given or 'all' -> process all images found in test_dir
%


% save current working directory
org_dir = pwd;

% get list of models
cd(model_dir);
models = dir('*.mat');
disp('Models found: ');
dir('*.mat');
newline;

% get list of test images
cd(org_dir);
cd(test_dir);
tests = dir('*.mat');
disp('Test images found: ');
dir('*.mat');
newline;

% parse arguments
if nargin < 9
  skip_done = true;
end
if nargin < 10 || strcmpi(subset, 'all')
  subset_ixs = 1:length(tests);             % by default all images in test dir are processed
  subset = 'all';
end


% process all pairs of model and test image
for m = models'
  m_base = m.name(1:(length(m.name)-4));
  disp(['Loading model ' m.name]);
  cd(org_dir);
  model = load([model_dir '/' m.name]);
  if strcmp(subset,'non-trn')
    if isfield(model.obj, 'trn_ixs')
      subset_ixs = setdiff(1:length(tests), model.obj.trn_ixs);  % only process non-trn images
      disp(['requested to process only non-training images; test subset = ' num2str(subset_ixs)]);
      newline;
    else
      error([mfilename ': requested to process only non-training images, but no .trn_ixs field found in model.']);
    end
  end
  %
  for t = tests(subset_ixs)'
    t_base = t.name(1:(length(t.name)-4));
    disp(['Detecting ' m_base ' in ' t_base]);
    disp('------------------------------------------');
  
    % should we process it ?
    doit = false;
    cd(org_dir);
    cd(process_dir);
    if ~skip_done
      doit = true;
    else
      process_fname = ['process_' m_base '_' t_base '.mat'];
      fid = fopen(process_fname);
      doit = (fid==-1);
      if fid>0
        fclose(fid);
        disp('process-file found -> case already processed or being processed -> skipping processing');
      end
    end
  
    % process
    if doit
      % write fake process* file to support parallalel processing
      trash=0; save(process_fname, 'trash');
      %
      % load data
      disp('Loading test datafile');
      cd(org_dir);
      test = load([test_dir '/' t.name]);
      %
      % compute PAS data on the fly if necessary
      if not(isfield(test.obj,'pas')) || not(isfield(test.obj,'ES'))
        cd(test_dir); % make it find _edges.tif
        test.obj = AddkASData(test.obj, false, true, 2);
      end
      %
      % detection
      disp('Detecting...');
      % recompute PASDall when matching a model from a different class than the one contained in the test image
      if isfield(det_pars,'recomp_pasdall') && det_pars.recomp_pasdall
          if isfield(test.obj,'PASDall')
            test.obj = rmfield(test.obj,'PASDall');
          end
      end
      H = ImagePASDATPSMatcher(model.obj, test.obj, HGA, det_pars, true);
      H = RecomputeNrgsAndScores(H, score_fct);
      % does all post-processing, including exporting figures and saving process file (last param == true)
      ConvertAllPASDATPS(model.obj, test.obj, H, test_dir, process_dir, draw_pars, remove_pars, true);
    end  % if doit

    newline;
  end  % loop over test images
end % loop over model images

% restore current working directory
cd(org_dir);

