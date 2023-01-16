%% analyseForagingModels 
% % inferential stats on experiment data
% poor = [subMean_poor, ones([length(subMean_poor), 1])];
% rich = [subMean_rich, 2*ones([length(subMean_rich), 1])];
%     y = array2table([poor; rich]);
%     tmp = stack(y, 1:3);
%     g1 = tmp.Var4;
%     g2 = tmp.Var1_Var2_Var3_Indicator;
%     [p, tbl, stats] = anovan(tmp.Var1_Var2_Var3,{g1, g2}, 'model', 'interaction');

