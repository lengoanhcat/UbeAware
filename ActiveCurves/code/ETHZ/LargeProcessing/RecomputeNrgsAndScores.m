function AH = RecomputeNrgsAndScores(AH, score_fct)

% Compute the cues-to-energy, and energy-to-score transformations
% according to score_fct(1) and score_fct(2) respectively.
%
% Every hypothesis AH(ix).H(hix) is scored by this function.
%


% compute energy from cues
switch score_fct(1).type

case 'weighed_sum'
  plain_nmp = strcmp(score_fct(1).note,'plain_nmp');
  if plain_nmp,  add_str = ' (plain NMP weight)';  else  add_str = '';  end
  disp(['Recompute energy with weigh sum, weights = [' num2str(score_fct(1).params,-1) ']' add_str]);
  AH = RecomputeNrg(AH, score_fct(1).params, plain_nmp);

case 'svm'
  AH = SVMScoreHyps(AH, score_fct(1).svm);

case 'tree'
  AH = DecTreeScoreHyps(AH, score_fct(1).T);

otherwise
 error([mfilename ': unknown energy function ' score_fct(1).type]);

end



% transform energy to scores
if length(score_fct) < 2
  error([mfilename ': you must provide an energy-to-score conversion function !']);
end


switch score_fct(2).type

case 'linear'
  nrg_roof = score_fct(2).params(2);  
  nrg_mult = score_fct(2).params(1);                              
  disp(['Recompute score as ' num2str(nrg_mult) '*nrg + ' num2str(nrg_roof)]);
  AH = LinNrg2Score(AH, nrg_mult, nrg_roof);   
  %AH = LinNrg2Score(AH, 1, 0);            % svm and decision tree scoring (-> .s = .nrg)

case 'gauss'
  stdr_dev = score_fct(2).params(1);
  disp(['Recompute score as exp(-(nrg^2)/(' num2str(stdr_dev) '^2))']);
  AH = GaussNrg2Score(AH, stdr_dev);       % score = exp(-(nrg^2)/(stdr_dev^2))

otherwise
  error([mfilename ': unknown energy-to-score conversion function ' score_fct(2).type]);

end

%AH = ProdScore(AH, 'basis', 4);            % bias score by initial Hough peak strength
