function AHC = ConvertAllPASDATPS(model, I, AH, test_dir, process_dir, draw_pars, remove_pars, save_process)

% Saves PASDATPS results,
% stored as sets of all hyps AH(ix).H(hix) found on all images I(ix),
% into standard process_{model_name}_{test_name}.mat file format containing all hyps H,
% their clustering C, and which clusters are retained in dets.
%
% WARNING: AH(ix).H must be sorted by .s (= score of hypothesis) !
% The higher the .s the better the hypothesis.
%

% save current working directory
org_dir = pwd;

% parse arguments
if nargin < 7
  save_process = false;
end
if nargin < 6
  draw_pars = false;
end

% process all pairs
proc_count = 0;
m_base = model.name;
for ix = 1:length(I)

  % flag info and basic shortcuts
  cd(org_dir); cd(process_dir);
  test = I(ix);
  t_base = test.name;
  proc_count = proc_count + 1;
  disp(['Test image ' num2str(proc_count) '/' num2str(length(I)) ' : ' t_base]); 
 
  % format hyps
  H = AH(ix).H;
  last_hix = 0;  
  for hix = 1:length(H)
    if H(hix).s >= 0
      last_hix = hix;
    end
    H(hix).m = [];                                        % no contour segments matched (is this field necessary?)
    if isfield(H(hix),'vx')                               % support Hough-only mode
      H(hix).bb = boundrect(H(hix).vx)';                  % hypothesis bounding-box
    elseif not(isfield(H(hix),'bb'))
      error([mfilename ': hypotheses must have either a .vx field, or a .bb field']);
    end
    H(hix).ms = sqrt(BBArea(H(hix).bb)/BBArea(model.bb)); % hypothesis scale
    H(hix).mT = mean(H(hix).bb')' - mean(model.bb')';     % hypothesis translation
    H(hix).mr = 0;                                        % no hypothesis rotation considered
  end

  % cluster hyps
  C = GeomGroupHyps(H(1:last_hix), model);                % includes summing score over hyps 
  for cix = 1:length(C)                                   % switch on for no hyp score accumulation
    C(cix).s = C(cix).simple_sc;
  end

  % filter hypothesis clusters
  [dets C] = RemoveContradictoryClusters(C, remove_pars);

  % automatic corr/wrong detection counter expects clusters sorted by their .s field
  % (this is not the case yet, as C is sorted by .simple_sc in GeomGroupHyps)
  [C new_order] = SortGroups(C, 's');
  [trash new_order_order] = sort(new_order);
  dets = new_order_order(dets);
  dets = sort(dets);

  % save process-file
  if save_process
    process_fname = ['process_' m_base '_' t_base '.mat'];
    disp(['Saving ' process_fname]);
    save(process_fname,'H','C','dets');
  end

  % draw and export figures with detections
  if not(islogical(draw_pars))
    cd(org_dir);
    if isempty(find(I(ix).ifname=='/',1)),  I(ix).ifname =  [test_dir  '/' I(ix).ifname];  end
    cd(process_dir);
    cids = ThresholdDetections(C(dets), 0.15);      % 0.15 -> for CVPR07 slides and paper figs
    if isempty(cids) && ~isempty(dets)
      cids = 1;                                     % keep at least one det (so draw every image)
    end
    %cids = ThresholdDetections(C(dets), 0.0);       % 0.04 -> for swan movies for CVPR07 poster (try to display everything)
    %cids = cids(1);                                 % only draw the 1 top-scored detection
    [trash fig_hs] = DrawTPSHyps(I(ix), H, false, false, draw_pars, C(dets(cids)));
    %ExportFigures(fig_hs, ['cool' m_base '_' t_base], 'jpg');        % use 'cool' to avoid overwrite (b/w for sup mat, 8.12.06)
    ExportFigures(fig_hs, ['forslides_' m_base '_' t_base], 'jpg');  % green for slides (6.2.07)
    %ExportFigures(fig_hs, ['forpaper_' m_base '_' t_base], 'jpg');   % for camready (29.4.07) and journal
    close(fig_hs);
  end

  % output
  AHC(ix).H = H;
  AHC(ix).C = C;
  AHC(ix).dets = dets;

  newline;

end  % loop over test images


% restore current working directory
cd(org_dir);
