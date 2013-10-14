function bb = ImageBB(I)

% bounding-box of the non-0 area of image I
%
% bb = [minx maxx miny maxy];
%

minx = find(sum(I,1),1);
miny = find(sum(I,2),1);
maxx = find(sum(I,1),1,'last');
maxy = find(sum(I,2),1,'last');

bb = [minx maxx miny maxy];
