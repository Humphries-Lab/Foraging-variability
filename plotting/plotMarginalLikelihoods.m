function plotMarginalLikelihoods(NLLEval, AlphaVector, BetaVector)

% plot parameter distributions - marginal likelihoods

numTiles = size(NLLEval, 2);
BlockTypes = {'rich', 'poor'};

for block = 1:2 % [rich poor]

    figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

    for n = 1:numTiles
        nexttile;
        SubjNLLEval = NLLEval{n}{block};
        lik = exp(-SubjNLLEval); % compute the likelihood, (rather than the log (negative LL))
        lik(isnan(lik)) = 0; % replace NANs with 0
        patch_tmp = sum(sum(lik, 3),2); % sum ml for alpha_patch across all beta and alpha_rho (MARGINAL likelihood)
        alpha_patch_marglik = patch_tmp/sum(patch_tmp);
        rho_tmp = sum(sum(lik, 3),1); % do the same for alpha_rho
        alpha_rho_marglik = rho_tmp/sum(rho_tmp);
        plot(AlphaVector, alpha_patch_marglik,AlphaVector, alpha_rho_marglik)
        xlabel('alpha', 'FontSize', 12)
        ylabel('p(data|model)', 'FontSize', 12)
        ylim([0 1])
        title(sprintf('Subject %d', n))

    end
    leg = legend({'Patch learning rate', 'Rho learning rate'}, 'FontSize', 14);
    leg.Layout.Tile = 'north';
    title(tl, sprintf('marginal likelihood distribution - learning rates - %s block', BlockTypes{block}), 'FontSize', 16)

    figure; tl = tiledlayout('flow', 'TileSpacing', 'Compact');

    for n = 1:numTiles
        nexttile;
        SubjNLLEval = NLLEval{n}{block};
        lik = exp(-SubjNLLEval); % compute the likelihood, (rather than the log (negative LL))
        lik(isnan(lik)) = 0; % replace NANs with 0
        beta_tmp = sum(squeeze(sum(lik, 1)),1); % do the same for betas
        beta_marglik = beta_tmp/sum(beta_tmp);
        plot(BetaVector, beta_marglik)
        xlabel('beta', 'FontSize', 12)
        ylim([0 1])
        ylabel('p(data|model)', 'FontSize', 12)
        title(sprintf('Subject %d', n))

    end
    leg = legend({'Softmax temperature' }, 'FontSize', 14);
    leg.Layout.Tile = 'north';
    title(tl, sprintf('marginal likelihood distribution - softmax temperature - %s block', BlockTypes{block}), 'FontSize', 16)
end

