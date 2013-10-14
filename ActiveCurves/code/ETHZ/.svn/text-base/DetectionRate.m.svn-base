function [det_rate, false_pos] = DetectionRate(S, C, t)

% detection rate and false positive rate
% as a function of the threshold t.
% S contans the scores and C the correct answer
% (1 for detected, 0 for not-detected)
%

% Flatten S and C to a simple list
S = reshape(S, 1, prod(size(S)));
C = reshape(C, 1, prod(size(C)));

% Compute real positive
D = (S >= t) .* C;

% Compute detection rate
det_rate = sum(D) / sum(C);

% Compute false positive
F = (S >= t) - D;

% Compute false positive rate
false_pos = sum(F) / sum(1-C);
