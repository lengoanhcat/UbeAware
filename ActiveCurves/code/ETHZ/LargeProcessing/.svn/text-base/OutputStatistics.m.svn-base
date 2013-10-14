function [xs, ys] = OutputStatistics(S, out_dir, col, thickness, curve_type)

% Plots statistics S.
%
% Also supports curve_type = 'acc/fppi' -> plot accuracy S.corr_acc as a function of the FPPI rate
%

% parse params
if nargin < 5
  curve_type = 'dr/fppi';
end
if nargin < 4
  thickness = 2;
end
if nargin < 3
  col = [0 0 1];
end
save_out = false;
if nargin < 2
  save_out = false;
elseif not(isempty(out_dir))
  save_out = true;
end

org_dir = pwd;
S = AccumulateStats(S);


% general info
disp(['Number of test images:          ' num2str(S.tot_test_imgs)]);
disp(['Number of positive images:      ' num2str(S.tot_pos)]);
disp(['Number of negative images:      ' num2str(S.tot_neg)]);
disp(['Number of instances:            ' num2str(S.tot_corr)]);
%
disp(['Number of correct detections:   ' num2str(length(S.corr_scores))]);
if isfield(S,'corr_overlaps');
  disp(['Average detection BB overlap:   ' num2str(mean(S.corr_overlaps))]);
end
disp(['Maximal detection rate:         ' num2str(length(S.corr_scores)/S.tot_corr)]);      
%disp(['Correct scores:                 ']); disp(num2str(sort(corr_scores)));
%disp(['Mean/std corr scores:           ' num2str([mean(corr_scores) std(corr_scores)])]);
%newline;
disp(['Number of false positives:      ' num2str(length(S.false_scores))]);
%disp(['False scores:                   ']); disp(num2str(sort(false_scores)));
%disp(['Mean/std false scores =         ' num2str([mean(false_scores) std(false_scores)])]);
%disp(['Hit <return> to continue']); keyboard;  % allow to play around with the data

if S.tot_pos == 0
  newline;
  disp('Test set contains no positive image. Make no plot.');
  return;
end


% generic statistics
[det_count false_count threshs] = ROCCount(S.corr_scores, S.false_scores);
threshs = [threshs(1) threshs];            % det_count(i) = number of corr dets with score > threshs(i) (strictly bigger than !)
det_count = [det_count(1) det_count];      % these 3 lines are for drawing the right most point of the plot
false_count = [10000 false_count];         % add rightmost point

% acc/fppi plot
if strcmp(curve_type, 'acc/fppi')
  acc = AverageAccuracy(S.corr_scores, S.corr_acc, threshs);
  newline;
  disp('Displaying ACCvFPPI plot');
  %figure;                                  % comment out to overlay curve over current figure
  xs = false_count/S.tot_test_imgs;
  ys = acc;
  plot(xs, ys, 'Color', col, 'LineWidth', thickness); hold on;         % default representation
  axis([0 1.5 0 max(ys)]);
  xlabel('False-positives per image');
  ylabel('Average accuracy');
  title('Accuracy vs FPPI plot');
end


% dr/fp plot
if strcmp(curve_type, 'dr/fp')
  newline;
  disp('Displaying DvFP plot');
  figure; plot(false_count, det_count, 'b', 'LineWidth', 1); hold on;  % draw continuous line
  plot(false_count, det_count, 'ok', 'LineWidth', 1);                  % highlight discrete points
  axis([0 S.tot_test_imgs 0 S.tot_corr]);
  xlabel('Number of false-positives');
  ylabel('Number of correct detections');
  disp('Exporting DvFP plot figure');
  if save_out
    cd(out_dir);
    %ExportFigureFix(gcf, 'dvfp.tif');
    disp('Saving DvFP data');
    save('stats.mat','S');
  end
end  % draw dr/fp plot ?


% draw histogram of scores (cool, see distribution, and can verify the correctness of the roc plot)
%newline;
%disp(['Displaying score histograms']);
%DrawHistograms(S.corr_scores, S.false_scores);


% DRvFPPI plot
if strcmp(curve_type, 'dr/fppi')
  newline;
  disp('Displaying DRvFPPI plot');
  %figure;  % comment out to overlay curve to current figure
  xs = false_count/S.tot_test_imgs;
  ys = det_count/S.tot_corr;
  plot(xs, ys, 'Color', col, 'LineWidth', thickness); hold on;  % default representation
  axis([0 1.5 0 1]);
  xlabel('False-positives per image');
  ylabel('Detection rate');
  title('Detection-rate vs FPPI plot');
  %
  % info about det-rates at a few specific fppi
  [trash ix] =  min(abs(xs - 0.1));
  disp(['Detection-rate at 0.1 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 0.2));
  disp(['Detection-rate at 0.2 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 0.3));
  disp(['Detection-rate at 0.3 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 1/3));
  disp(['Detection-rate at 1/3 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 0.35));
  disp(['Detection-rate at 0.35 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 0.4));
  disp(['Detection-rate at 0.4 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  [trash ix] =  min(abs(xs - 0.5));
  disp(['Detection-rate at 0.5 FPPI = ' num2str(ys(ix),3) ' (thresh = ' num2str(threshs(ix)) ')']);
  %
  if save_out
    disp('Exporting DRvFPPI plot figure');
    ExportFigureFix(gcf, 'drvfppi.tif');
   end
end % draw dr/fppi plot ?


% precision/recall plot
if strcmp(curve_type, 'prec/recall')
  newline;
  disp('Displaying Precision/Recall plot');
  %[pr_eer th] = PrecRecallPlot(S.corr_scores, S.false_scores, S.tot_corr, true, true, col); % first logical -> new fig, second logical -> draw eer line
  pr_eer = PrecRecallPlot(S.corr_scores, S.false_scores, S.tot_corr, false, true, col); % don't open new fig -> overlay existing one
  %
  disp(['Recall at equal error rate = ' num2str(pr_eer,3) ' (thresh = ' num2str(th) ')']);
  %
  if save_out
    disp('Exporting prec-recall plot figure');
    ExportFigureFix(gcf, 'prec-recall.tif');
  end
end % draw prec/recall plot ?


% ROC plot
if strcmp(curve_type, 'roc')
  newline;
  if S.tot_neg > 1 && S.tot_pos > 1
    disp('Consider the task as a binary classification problem on presence/absence of the object class');
    disp('Displaying ROC plot');
    roc_eer = ROCPlot([S.pos_scores S.neg_scores], [ones(1,S.tot_pos) zeros(1,S.tot_neg)], true, true);
    %
    disp(['Detection rate at equal error rate = ' num2str(roc_eer,3)]);
    %
    if save_out
      cd(org_dir);
      disp('Exporting ROC plot figure');
      ExportFigureFix(gcf, 'roc.tif');
    end
  else
    disp('The test set is not composed of a mixture of positive and negative images -> cannot plot classification ROC');
  end
end % display roc plot ?


% restore original directory
cd(org_dir);
