function [allfilter, allsymbol] = makefilter(scale, norient) 
% make Gabor filters at fixed scale with norient orientations
% "allfilter" contains all the Gabor sine and cosine pairs
% "allsymbol" contains bars for display purpose

allorient  = (0:norient-1) * 180/norient; 
allfilter = cell(1, norient); 
allsymbol = cell(1, norient); 
for o = 1 : norient
    [allfilter{1, o}, allsymbol{1, o}] = gaborfilter(scale, allorient(o),'single');
end
    



    

