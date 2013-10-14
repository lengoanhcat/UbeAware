function H = MakeGTOBBShapes(GTO, gto_type, gto_ixs)

% Generate a bounding-box as a shape H(ix) for each true shape GTO(:,:,gto_ixs)
%
%

GTO = SplitAndThinGTO(GTO, gto_type);

ix = 0;
for kix = gto_ixs
  %
  % determine bb
  bb = ImageBB(GTO(:,:,kix));   % bb = [minx maxx miny maxy]
  %figure; hold on;
  %imshow(GTO(:,:,kix));
  %DrawBB(bb, [1 0 0]);
  %
  % produce BB's 'shape'
  pts = [];
  minx = bb(1); maxx = bb(2);
  miny = bb(3); maxy = bb(4);
  pts = [pts [minx:maxx; ones(1,maxx-minx+1)*miny]];
  pts = [pts [minx:maxx; ones(1,maxx-minx+1)*maxy]];
  pts = [pts [ones(1,maxy-miny+1)*minx; miny:maxy]];
  pts = [pts [ones(1,maxy-miny+1)*maxx; miny:maxy]];
  %plot(pts(1,:), pts(2,:), 'go');
  %
  % format it as a hyp H
  ix = ix + 1;
  H(ix).vx = pts';
  H(ix).vy = pts';
end
