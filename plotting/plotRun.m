function plotRun(model, ExampleRun)
% figure - example run

SetUpEnviron

figure
plot(1:Env.BlockTime, ExampleRun.EstimatedPatchRR, 'LineWidth', 1) % plot patch reward rate
hold on
plot(1:Env.BlockTime, ExampleRun.Rho, 'LineWidth', 1) % plot estimated average RR
legend('Patch RR', '\rho: Estimated average RR', 'FontSize', 16, 'FontName', 'Helvetica')
xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Reward rates','FontSize', 16, 'FontName', 'Helvetica');

print(sprintf('~/Dropbox/foraging/outputs/M%d/M%d_example_run', model, model), '-dpng')