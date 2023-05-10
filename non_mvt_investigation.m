%% plots to show what results we have for MVT (non) model 

clear
close all

%load ~/Dropbox/foraging/outputs/M2/fitting_results_separate_230504.mat %
%for reference 
load ~/Dropbox/foraging/outputs/M2/fitting_results_separate.mat % latest one

% plot distribution of alpha patch
figure; 
boxchart(squeeze(minNLLFitParams(:,1,:)))
title('Best fit parameter values for alpha patchRR', 'FontSize', 18)
ylim([0,1])
xlabel('Environment', 'FontSize', 14); xticklabels({'Rich', 'Poor'});
ylabel('Alpha patch RR value', 'FontSize', 14)

% plot distribution of alpha rho
figure; 
boxchart(squeeze(minNLLFitParams(:,1,:)))
title('Best fit parameter values for alpha rho', 'FontSize', 18)
ylim([0,1])
xlabel('Environment', 'FontSize', 14); xticklabels({'Rich', 'Poor'});
ylabel('Alpha rho value', 'FontSize', 14)

% plot distribution of beta
figure; 
boxchart(squeeze(minNLLFitParams(:,2,:)))
title('Best fit parameter values for beta', 'FontSize', 18)
xlabel('Environment', 'FontSize', 14); xticklabels({'Rich', 'Poor'});
ylabel('Beta value', 'FontSize', 14)

