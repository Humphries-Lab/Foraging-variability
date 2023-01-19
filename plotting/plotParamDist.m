function plotParamDist(FitParams)

% parameter distribution for each subject
ParamNames = {'AlphaPatch-Q', 'AlphaRho', 'InverseTemp'};
numTiles = length(FitParams);
BlockTypes = {'rich', 'poor'};

for block = 1:2 % [rich poor]

    for iParam = 1:3
        figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');
        title(tl, sprintf('Distribution of parameter fits - %s - %s',ParamNames{iParam}, BlockTypes{block}), 'FontSize', 16);

        for n = 1:numTiles
            nexttile;
            histogram(FitParams{n}{block}(:,iParam)); % alpha patch
            title(sprintf('Subject %d', n))
            xlabel(sprintf('%s - Fits',ParamNames{iParam}), 'FontSize', 12)
        end
    end
end