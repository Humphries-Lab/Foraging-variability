function plotBestParams(minNLLFitParams)
% plot the best parameters for each person
m = squeeze(median(minNLLFitParams));

% parameter distribution for each subject
%ParamNames = {'AlphaPatch-Q', 'AlphaRho', 'InverseTemp'};
ParamNames = {'AlphaRho', 'InverseTemp'};
BlockTypes = {'rich', 'poor'};

for block = 1:2 % [rich poor]

    figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
    title(tl, sprintf('Best fitting parameters across subjects - %s', BlockTypes{block}), 'FontSize', 16);

    nexttile
    plot(minNLLFitParams(:,1,block),minNLLFitParams(:,2,block),'o','MarkerSize',10); hold on
    plot(m(1,block),m(2,block),'+','MarkerSize',15)
    xlabel(ParamNames{1},'FontSize', 12);
    ylabel(ParamNames{2},'FontSize', 12);

%     nexttile
%     plot(minNLLFitParams(:,1,block),minNLLFitParams(:,3,block),'o','MarkerSize',10); hold on
%     plot(m(1,block),m(3,block),'+','MarkerSize',15)
%     xlabel(ParamNames{1},'FontSize', 12);
%     ylabel(ParamNames{3},'FontSize', 12);
% 
%     nexttile
%     plot(minNLLFitParams(:,2,block),minNLLFitParams(:,3,block),'o','MarkerSize',10); hold on
%     plot(m(2,block),m(3,block),'+','MarkerSize',15)
%     xlabel(ParamNames{2}, 'FontSize', 12);
%     ylabel(ParamNames{3},'FontSize', 12);

end

leg = legend({'Subject fit', 'Group Median'}, 'FontSize', 14);
leg.Layout.Tile = 'north';
