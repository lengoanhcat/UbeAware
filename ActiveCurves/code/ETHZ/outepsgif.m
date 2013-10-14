function outepsgif(pic, outname); 

[y x c] = size(pic); 
figure('Units','Pixels','Resize','off',...
   'Position',[100 100 x y],'PaperUnit','points',...
   'PaperPosition',[0 0 x y]);
axes('position',[0 0 1 1]);
imshow(pic, []);
axis off
saveas(gcf, outname, 'eps');

towrite = double(pic);  
towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
imwrite(towrite, [outname '.gif'],'gif');