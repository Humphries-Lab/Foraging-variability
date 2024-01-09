
figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(1,mean(minNLLFitParams_rich.beta),'FaceColor',color.rich, 'EdgeColor',color.rich); hold on 
bar(2,mean(minNLLFitParams_poor.beta),'FaceColor',color.poor, 'EdgeColor',color.poor); 

% plot individual points
scatter(ones(size(minNLLFitParams_rich.beta)), minNLLFitParams_rich.beta,'o','MarkerFaceColor', marker.rich,'MarkerEdgeColor', marker.rich,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 
scatter(2*ones(size(minNLLFitParams_poor.beta)), minNLLFitParams_poor.beta,'o','MarkerFaceColor', marker.poor,'MarkerEdgeColor', marker.poor,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 

% plot error bars 
betaSEM_rich = std(minNLLFitParams_rich.beta)./sqrt(numel(minNLLFitParams_rich.beta));
betaSEM_poor = std(minNLLFitParams_poor.beta)./sqrt(numel(minNLLFitParams_poor.beta));

errorbar(1,mean(minNLLFitParams_rich.beta),betaSEM_rich,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 
errorbar(2,mean(minNLLFitParams_poor.beta),betaSEM_poor,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 

set(gca,'XTick',1:2, 'XTickLabel', {'Rich', 'Poor'}, 'YScale', 'log', 'YLim', [0,1])
ylabel('Fit beta (higher = exploit)')
title('Model 26')

% statistics
[pearson_coeff, p_val]  = corr([minNLLFitParams_rich.beta; minNLLFitParams_poor.beta], [meanLT_rich'; meanLT_poor'], "type","Spearman")
[pearson_coeff, p_val]  = corr([minNLLFitParams_rich.bias; minNLLFitParams_poor.bias], [meanLT_rich'; meanLT_poor'], "type","Spearman")

[p_beta, h, t] = signrank(minNLLFitParams_rich.beta, minNLLFitParams_poor.beta)
[p_beta, h, t] = signrank(minNLLFitParams_rich.alphaRho, minNLLFitParams_poor.alphaRho)
[p_beta, h, t] = signrank(minNLLFitParams_rich.bias, minNLLFitParams_poor.bias)


figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(1,mean(minNLLFitParams_rich.alphaRho),'FaceColor',color.rich, 'EdgeColor',color.rich); hold on 
bar(2,mean(minNLLFitParams_poor.alphaRho),'FaceColor',color.poor, 'EdgeColor',color.poor); 

% plot individual points
scatter(ones(size(minNLLFitParams_rich.alphaRho)), minNLLFitParams_rich.alphaRho,'o','MarkerFaceColor', marker.rich,'MarkerEdgeColor', marker.rich,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 
scatter(2*ones(size(minNLLFitParams_poor.alphaRho)), minNLLFitParams_poor.alphaRho,'o','MarkerFaceColor', marker.poor,'MarkerEdgeColor', marker.poor,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 

% plot error bars 
alphaRhoSEM_rich = std(minNLLFitParams_rich.alphaRho)./sqrt(numel(minNLLFitParams_rich.alphaRho));
alphaRhoSEM_poor = std(minNLLFitParams_poor.alphaRho)./sqrt(numel(minNLLFitParams_poor.alphaRho));

errorbar(1,mean(minNLLFitParams_rich.alphaRho),alphaRhoSEM_rich,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 
errorbar(2,mean(minNLLFitParams_poor.alphaRho),alphaRhoSEM_poor,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 

set(gca,'XTick',1:2, 'XTickLabel', {'Rich', 'Poor'})
ylabel('Fit alphaRho')
title('Model 5')




figure('Units', 'centimeters', 'PaperPositionMode', 'auto' ,'Position',figsize.square);
bar(1,mean(minNLLFitParams_rich.bias),'FaceColor',color.rich, 'EdgeColor',color.rich); hold on 
bar(2,mean(minNLLFitParams_poor.bias),'FaceColor',color.poor, 'EdgeColor',color.poor); 

% plot individual points
scatter(ones(size(minNLLFitParams_rich.bias)), minNLLFitParams_rich.bias,'o','MarkerFaceColor', marker.rich,'MarkerEdgeColor', marker.rich,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 
scatter(2*ones(size(minNLLFitParams_poor.bias)), minNLLFitParams_poor.bias,'o','MarkerFaceColor', marker.poor,'MarkerEdgeColor', marker.poor,'MarkerEdgeAlpha', M_alpha,'MarkerFaceAlpha', M_alpha); 

% plot error bars 
biasSEM_rich = std(minNLLFitParams_rich.bias)./sqrt(numel(minNLLFitParams_rich.bias));
biasSEM_poor = std(minNLLFitParams_poor.bias)./sqrt(numel(minNLLFitParams_poor.bias));

errorbar(1,mean(minNLLFitParams_rich.bias),biasSEM_rich,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 
errorbar(2,mean(minNLLFitParams_poor.bias),biasSEM_poor,'LineStyle', 'none','LineWidth',widths.plot,'Color',[0.1 0.1 0.1],'CapSize', 0); 

set(gca,'XTick',1:2, 'XTickLabel', {'Rich', 'Poor'})
ylabel('Fit bias (nmore negative = overharvesting')
title('Model 26')