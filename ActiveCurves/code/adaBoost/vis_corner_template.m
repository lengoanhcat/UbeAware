% visulize selected corners
%% decode corner position and orientations
corner_index = [];
for i_feature = 1:length(model)
    dim = model(i_feature).dimension;
    dim = p_map(dim)-p1+1;
    index = find(m3VecBound(:,3)<dim,1,'last');
    orients = m3VecBound(index,1:2);
    residual_dim = dim - m3VecBound(index,3);
    [posX posY] = ind2sub([sx sy],residual_dim);
    corner_index = [corner_index;[orients,posX,posY]];
end
% draw corners
ORI_I = corner_index(:,1);
ORI_J = corner_index(:,2);
X = corner_index(:,3);
Y = corner_index(:,4);
 figure;
 
    imshow(ones(sx,sy));
    hold on
    line([1, 1,sy,sy,1],[1,sx,sx,1,1],'Color','k')
    for i = 1:length(model)
        len = 10;
        alpha = (ORI_I(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'r.-','linewidth',1);
        
        alpha = (ORI_J(i)-1)*pi/nOrient;
        plot([Y(i) Y(i)+round(len*cos(alpha))], [X(i) X(i)-round(len*sin(alpha))],'g.-','linewidth',1);
        set(gcf,'PaperPositionMode','auto');
        imgName= sprintf('../results/CN_template_Fold_%d_feature_%d.png',iFold,i);
        saveas(gcf,imgName);
    end

  