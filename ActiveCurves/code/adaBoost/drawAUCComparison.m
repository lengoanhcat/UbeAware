% code to draw the AUC comparison between three feature groups
fileName = [ outPath '/Fold_5_result.mat'];
load(fileName);
AB_auc_table = aucTable;


fileName = [ outPath '/ARC_Fold_5_result.mat'];
load(fileName);
arc_auc_table = aucTable;

fileName = [ outPath '/AP_Fold_5_result.mat'];
load(fileName);
corner_auc_table = aucTable;

figure;
hold on;
nFeature = size(AB_auc_table,2);
errorbar(1:nFeature,mean(AB_auc_table),std(AB_auc_table),'b','linewidth',2)

nFeature = size(arc_auc_table,2);
errorbar(1:nFeature,mean(arc_auc_table),std(arc_auc_table),'g','linewidth',2);

nFeature = size(corner_auc_table,2);
errorbar(1:nFeature,mean(corner_auc_table),std(corner_auc_table),'r','linewidth',2);

axis([0 nFeature 0.5 1])
legend('active basis','active arcs','active arcs + corners','Location','SouthEast')
xlabel('number of featuers')
ylabel('AUC')

figure;
hold on;
nFeature = size(AB_auc_table,2);
plot(1:nFeature,mean(AB_auc_table),'bo-','linewidth',2)

nFeature = size(arc_auc_table,2);
plot(1:nFeature,mean(arc_auc_table),'gx-','linewidth',2);

nFeature = size(corner_auc_table,2);
plot(1:nFeature,mean(corner_auc_table),'r+-','linewidth',2);

axis([0 nFeature 0.5 1])
legend('active basis','active arcs','active arc + corners','Location','SouthEast')
xlabel('number of featuers')
ylabel('AUC')

fileName = [ outPath '/auc_figure'];
saveas(gcf,[fileName '.fig']);
saveas(gcf,[fileName '.png']);
