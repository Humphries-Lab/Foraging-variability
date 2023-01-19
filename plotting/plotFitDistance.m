function plotFitDistance(FitParams, StartParams, nStarts)

figure
scatter3(FitParams(:,1),FitParams(:,2),FitParams(:,3), 50, FitParams(:,3), 'filled');
hold on;
scatter3(StartParams(:,1),StartParams(:,2),StartParams(:,3), 50, StartParams(:,3), "r", 'filled');
colormap("winter");
colorbar;
xlabel('alpha patch', 'FontSize', 12)
ylabel('alpha rho', 'FontSize', 12)
zlabel('beta', 'FontSize', 12)
% plot lines between the start and fit points
for ii = 1:nStarts
    v1 = [StartParams(ii,1),StartParams(ii,2),StartParams(ii,3)];
    v2 = [FitParams(ii,1),FitParams(ii,2),FitParams(ii,3)];
    v = [v2;v1];
    plot3(v(:,1), v(:,2),v(:,3),'r', 'LineStyle', '--')
end

title('Distance from start to fit parameters', 'FontSize', 16)
legend({'Fit parameters', 'Start parameters'}, 'FontSize', 14);
end