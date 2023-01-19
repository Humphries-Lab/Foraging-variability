function plotRun(model)
% figure - example run

SetUpEnviron
load(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model))

figure
plot(1:Env.BlockTime, out.PatchRR, 'LineWidth', 1) % plot patch reward rate
hold on
plot(1:Env.BlockTime, out.Rho*57.5, 'LineWidth', 1) % plot estimated average RR
legend('Patch RR', '\rho: Estimated average RR', 'FontSize', 16, 'FontName', 'Helvetica')
xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Reward rates','FontSize', 16, 'FontName', 'Helvetica');

print(sprintf('~/Dropbox/foraging/outputs/M%d/M%d_example_run', model, model), '-dpng')