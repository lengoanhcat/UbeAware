function [GTO, gto_type] = LoadGTOutlines(test_dir, t_base, class_name)

% Vitto format
gto_fname = [test_dir '/' t_base '_' class_name '_outlines.pgm'];
if exist(gto_fname,'file')
  disp(['Loading ground-truth outlines ' gto_fname]);
  GTO = imread(gto_fname);
  gto_type = 'vitto';
  return;
end

% Marcin format
gto_fnames = dir([test_dir '/' t_base '.mask*.png']);
if not(isempty(gto_fnames))
  k = 0;
  for gto_fname = gto_fnames'
    k = k + 1;
    disp(['Loading ground-truth outlines ' gto_fname.name]);
    GTO(:,:,k) = imread([test_dir '/' gto_fname.name]);
  end        
  gto_type = 'marcin';
  return
end

% no gto found
GTO = [];
gto_type = false;
