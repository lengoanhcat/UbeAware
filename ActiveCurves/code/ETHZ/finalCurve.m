addpath('./LargeProcessing');
colors(1,:) = [1 1 0];           % yellow
colors(2,:) = [0 1 0];           % green
colors(3,:) = [1 0 0];           % red
colors(4,:) = [0 0 1];           % blue
colors(5,:) = [0 1 1];           % cyan
colors(6,:) = [0.6 0.6 0.6];     % grey
colors(7,:) = [1 0 1];           % purple
colors(8,:) = [1 0.5 0];         % orange
colors(9,:) = [0 0.2 0.6];       % dark blue
colors(10,:) = [0 0.3 0];        % dark green
colors(11,:) = [1 1 1];          % white
colors(12,:) = [1 0.2 0.6];      % brown
colors(13,:) = [0 0 0];          % black
colors(14,:) = [0.3 0 0.7];      % blue-purple
colors(15,:) = [0.5 0.5 0];      % dark yellow
colors(16,:) = [1 0.6 0.8];      % bright brown

%%%%%

result_path = '../../results/ETHZ';
dirs = dir([ result_path '/result*']);
clname = {'Applelogos', 'Bottles', 'Giraffes', 'Mugs', 'Swans', 'FrontalFaces2'};   % classes considered

figure; hold on;
Nmodels  = length(dirs);
for i_cl = 1:6
    clear curves;
    subplot(1,6,i_cl);
    hold on;
    title(clname{i_cl});
    % curve for each run
    for model_ix = 1:Nmodels
        load([result_path '/' dirs(model_ix).name '/' clname{i_cl} '_SV2.mat']);
        [xs ys] = OutputStatistics(S, [], colors(model_ix,:), 1, 'dr/fppi');
        curves{model_ix} = [xs; ys];
    end
    % curve fore all the 5 runs
    if Nmodels > 1
        [mc stdvs] = FuseCurves(curves, 1.5, false);   % 1.5 -> max FPPI val; last param -> display
        plot(mc(1,:),mc(2,:),'color',colors(16,:),'lineWidth',3);
        target_gcf=hgload(['./real-images/' clname{i_cl} '.fig']);
        plot(mc(1,:),mc(2,:),'color',colors(3,:),'lineWidth',2);
        axis([0 1.5 0 1]);
        saveas(gcf,[ result_path '/' clname{i_cl} '.fig']);
        saveas(gcf,[ result_path '/' clname{i_cl} '.pdf']);
        close;
        newline;
        disp(['Mean accuracy / stdv at 0.4 FPPI: ' num2str([mc(2,9) stdvs(9)])]);
        newline;
    end
end