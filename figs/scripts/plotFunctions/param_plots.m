

meanLT_rich = mean(squeeze(subject_leave_times.mean(1,:,:)));
meanLT_poor = mean(squeeze(subject_leave_times.mean(2,:,:)));

figure
scatter(minNLLFitParams.bias_rich, meanLT_rich,  70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w'); hold on
scatter(minNLLFitParams.bias_poor, meanLT_poor,  70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w'); hold on
%set(gca,'XLim', [0 2],'XScale', 'log','XTick', [10^-2, 10^-1, 10^0, 10^1])
xlabel('Fit bias (higher = exploit)')
ylabel('Patch leaving time (s)')
legend({'rich', 'poor'})
set(findall(gcf,'-property','FontSize'),'FontSize',18)

tmp_set1 = brewermap(9,'Set1');
color.rich = tmp_set1(2,:); % green
color.poor = tmp_set1(3,:); % purple


figure
scatter(minNLLFitParams.bias_rich, minNLLFitParams.beta_rich, 70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w'); hold on
scatter(minNLLFitParams.bias_poor, minNLLFitParams.beta_poor, 70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w')

plot([minNLLFitParams.bias_rich(:)';minNLLFitParams.bias_poor(:)'], [minNLLFitParams.beta_rich(:)';minNLLFitParams.beta_poor(:)'], 'k-')
legend({'rich', 'poor'})
xlabel('Bias')
ylabel('Beta')
set(findall(gcf,'-property','FontSize'),'FontSize',18)

%% environment differences in beta fits 

 tmp_set1 = brewermap(9,'Set1');
 color.rich = tmp_set1(2,:); % green
  color.poor = tmp_set1(3,:); % purple

figure
scatter(1, minNLLFitParams.beta_rich, 70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w'); hold on
scatter(2, minNLLFitParams.beta_poor, 70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w')

ylabel('Beta')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'}, 'YScale', 'log')
errorbar(1, mean(minNLLFitParams.beta_rich), std(minNLLFitParams.beta_rich)./sqrt(numel(minNLLFitParams.beta_rich)), "both", "CapSize", 0, "LineWidth", 3, "Color", 'k')
errorbar(2, mean(minNLLFitParams.beta_poor), std(minNLLFitParams.beta_poor)./sqrt(numel(minNLLFitParams.beta_rich)), "both", "CapSize", 0, "LineWidth", 3, "Color", 'k')

%% bias fits - no environment
figure
scatter(1, minNLLFitParams.bias, 'k', 'filled'); hold on

ylabel('Bias')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(gca,'XLim', [0 2],'XTick', [1], 'XTickLabel', {'Bias fit'})
errorbar(1, mean(minNLLFitParams.bias), std(minNLLFitParams.bias)./sqrt(numel(minNLLFitParams.bias)), "both", "CapSize", 0, "LineWidth", 3, "Color", 'k')


%% bias fits - environment
figure
scatter(1, minNLLFitParams.bias_rich, 70, 'MarkerFaceColor',color.rich, 'MarkerEdgeColor','w'); hold on
scatter(2, minNLLFitParams.bias_poor, 70, 'MarkerFaceColor',color.poor, 'MarkerEdgeColor','w')

ylabel('bias')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(gca,'XLim', [0 3],'XTick', [1:2], 'XTickLabel', {'Rich', 'Poor'})
errorbar(1, mean(minNLLFitParams.bias_rich), std(minNLLFitParams.bias_rich)./sqrt(numel(minNLLFitParams.bias_rich)), "both", "CapSize", 0, "LineWidth", 3, "Color", 'k')
errorbar(2, mean(minNLLFitParams.bias_poor), std(minNLLFitParams.bias_poor)./sqrt(numel(minNLLFitParams.bias_rich)), "both", "CapSize", 0, "LineWidth", 3, "Color", 'k')
