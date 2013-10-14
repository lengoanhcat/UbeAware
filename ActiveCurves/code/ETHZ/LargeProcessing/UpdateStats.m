function [S, missed, corr_dets] = UpdateStats(S, gt_fname, C, dets, keep, pascal)

% Update statistics S with new evidence C.
%
% If 'keep' given: only keep the top 'keep' dets (if available)
%
% Returns also corr_dets, which is an index into C: C(corr_dets) = correct detections
% (hence, corr_dets is a subset of dets)
%
% parse arguments
if nargin < 6
  pascal = false;
end
if nargin >= 5 && not(strcmp(keep,'all'))
  keep = min(keep,length(dets));
  dets = dets(1:keep);
end

if pascal
  disp('Using PASCAL criterion');
else
  disp('Using Ferrari criterion');
end


disp(['Attempting to load groundtruth file ' gt_fname]);
gt_found = exist(gt_fname, 'file');
if not(gt_found == 2)
  % try without s
  dl = find(gt_fname == '.'); % location of '.'
  if gt_fname(dl-1) == 's'
    disp('Removing ''s'' from groundtruth filename and trying again');
    gt_fname = [gt_fname(1:(dl-2)) gt_fname(dl:end)];
    gt_found = exist(gt_fname, 'file');
  end
end


if gt_found

  % a positive image
  disp('Positive image');
  S.tot_pos = S.tot_pos+1;
  if ~isempty(dets)
    S.pos_scores = [S.pos_scores C(dets(1)).s];
  else
    S.pos_scores = [S.pos_scores -inf];
  end
  gtBB = load(gt_fname);
  [corr_dets false_pos cs fs corr_overlaps] = CorrectDetections(C(dets), gtBB, pascal);       
  S.corr_scores = [S.corr_scores cs];
  S.false_scores = [S.false_scores fs];
  S.tot_corr = S.tot_corr+size(gtBB,1);
  S.corr_overlaps = [S.corr_overlaps corr_overlaps];
  missed = size(gtBB,1)-length(cs);
  disp([num2str(length(corr_dets)) ' correct detections : ' num2str(cs)]);
  disp(['overlaps             : ' num2str(corr_overlaps)]);
  disp([num2str(missed) ' missed detections']);
  disp([num2str(length(false_pos)) ' false-positives    : ' num2str(fs)]);
  disp(['score of image       : ' num2str(S.pos_scores(end))]);

else   % no groundtruth file -> no instances of the model class in this test image

  % a negative image
  disp('Negative image');
  S.tot_neg = S.tot_neg+1;
  corr_dets = [];
  disp('No groundtruth file found -> assume no instances in this image');
  if ~isempty(dets)
    S.neg_scores = [S.neg_scores C(dets(1)).s];
  else
    S.neg_scores = [S.neg_scores -inf];
  end
  fs = [];
  for j = dets
    fs = [fs C(j).s];
  end
  S.false_scores = [S.false_scores fs];
  disp([num2str(length(dets)) ' false-positives, ' num2str(fs)]);
  disp(['score of image ' num2str(S.neg_scores(end))]);
  missed = 0;

end % ground-truth file found

% set outputs
S.tot_test_imgs = S.tot_test_imgs+1;
corr_dets = dets(corr_dets);
