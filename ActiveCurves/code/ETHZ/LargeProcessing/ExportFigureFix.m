function ExportFigureFix(fig, fname, resize, crop)

% Exports figure fig to file fname.
% Fixes the output to have no border and correct resolution.
%
% resize = true if want to resize to true image size (exporting makes always 1200x901 images)
% crop = false if no additional crop; otherwise crop = [x1 y1 width height]
%
% File format is determined by last 3 chars of fname
%

if nargin < 3
  resize = true;
end

if nargin < 4
  crop = false;
end

% get real width and height of the figure
p = get(fig,'Position');
w = p(3); h = p(4);

% export figure;
% always has resolution 1200x901, with the real figure
% covering the whole height, and a white border proportional on the width.
% !!!! CAN USE -r0 OPTION TO KEEP SCREEN RESOLUTION !! YEEEAH!
if all(fname((end-2):end) == 'tif') || all(fname((end-3):end) == 'tiff')
  print(fig, '-dtiff', fname);
elseif all(fname((end-2):end) == 'jpg')
  print(fig, '-djpeg', fname);
elseif all(fname((end-2):end) == 'png')
  print(fig,'-dpng',fname);
else
  error(['ExportFigureFix: ' fname ' file type not supported']);
end

% fix border and resolution
f = 901/h;       % scale factor
x1 = (1200-w*f)/2;
x2 = 1200-x1;

% remove border
system(['convert -crop ' num2str(x2-x1) 'x901' '+' num2str(x1) '+0 "' fname '" "' fname '"']); 

if resize
  system(['convert -resize ' num2str(w) 'x' num2str(h) ' "' fname '" "' fname '"']);
end

if length(crop) > 1
  if not(resize)
    crop = crop * f;
  end
  system(['convert -crop ' num2str(crop(3)) 'x' num2str(crop(4)) '+' num2str(crop(1)) '+' num2str(crop(2)) ' "' fname '" "' fname '"']); 
end
