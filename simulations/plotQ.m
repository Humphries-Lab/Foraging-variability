function plotQ(model)

SetUpEnviron
load(sprintf('~/Dropbox/foraging/outputs/M%d/SimulationResults.mat', model))

figure
plot(1:Env.BlockTime, out.Q(:,2), 'LineWidth', 1) % plot QStay
hold on
plot(1:Env.BlockTime, out.Q(:,1), 'LineWidth', 1) % plot leave
legend('QStay', 'QLeave', 'FontSize', 16, 'FontName', 'Helvetica')
xlabel('Time (s)','FontSize', 16, 'FontName', 'Helvetica');
ylabel('Estimated value','FontSize', 16, 'FontName', 'Helvetica');

print(sprintf('~/Dropbox/foraging/outputs/M%d/M%d_example_run', model, model), '-dpng')