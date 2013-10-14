function S = InitStats 

S.corr_scores = [];           % score of every single correct detection (BB)
S.corr_acc = [];              % shape matching accuracy of corr dets
S.gt_acc = [];                % shape accuracy of ground-truth BB
S.corr_ids = [];              % [img_id; gt_id]: ground-truth outline associated to each correct detection
S.false_scores = [];          % same for false-pos
S.corr_overlaps = [];         % percentage area overlaps with ground-truth BB (mean over two dirs)
S.pos_scores = [];            % score of each positive image (= score of the highest scored dectection within it)
S.neg_scores = [];            % same for neg images
S.tot_corr = 0;               % number of instances of the class (>= tot_pos)
S.tot_test_imgs = 0;          % number of test images
S.tot_pos = 0;                % number of images containing >=1 instances of the class
S.tot_neg = 0;                % number of images not containing any instance of the class
