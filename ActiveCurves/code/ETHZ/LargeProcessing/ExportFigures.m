function ExportFigures(fig_hs, header, format)

% export figures fig_hs to pictures in the given format.
% Filenames are:
%   header_index.format   if   more than 1 fig, or
%   header.format         if   just one fig
%

ix = 0;
for h = fig_hs

  ix = ix+1;       % use ix and not h as filename extension, as fig_hs is not guaranteed to start with 1
  if length(fig_hs) > 1
    ExportFigureFix(h, [header '_' num2str(ix) '.' format]);
  else
    ExportFigureFix(h, [header '.' format]);
  end

end
