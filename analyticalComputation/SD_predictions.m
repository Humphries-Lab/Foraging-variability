%% analysis of SDs of leaving times for Le Heron data 
% Emma Scholey 4 May 2023
clear
close all;
addpath(genpath('~/Dropbox/foraging/code'))

blockNames = {'rich', 'poor'};

%% analysis of real subject data 
% get their data
subj = [1:23 25:40];

% initialise
numSubjects = size(subj,2);
numBlocks = size(blockNames, 2);

%subj = 1:4;
load ~/Dropbox/foraging/raw_data/summary/young_variables/t_young.mat

%% per-patch, per-environment SD 
% test model predictions 

SD_LT = zeros([numSubjects, 3, 2]); % each row is subject SD, each column patch type, 3d is block type
for iSubj = 1:numSubjects % for each subject
    for iPatch = 1:3
        for iBlock = 1:2
            SD_LT(iSubj,iPatch,iBlock) = std(t_young.leaveT(t_young.subj == iSubj & t_young.env == iBlock & t_young.patch == iPatch));
        end
    end
end

% between environments 

real_change_in_SD_between_environments = mean(SD_LT(:,:,2) - SD_LT(:,:,1)); % poor - rich
% our results agree in that there is higher SD in poor environment compared
% to rich, but differences in real SDs are much smaller than analytical SDs
% and don't see increase in SD differences as patch yield increases
% e.g. (0.5, 0.03 0.22) compared to (4.5, 7, 9)  

% within environments 

real_change_in_SD_within_rich = mean(SD_LT(:,:,1)) - mean(SD_LT(:,:,1))';

real_change_in_SD_within_poor = mean(SD_LT(:,:,2)) - mean(SD_LT(:,:,2))';

% difference between low and high yield patches in rich environment is
% 0.556, and in poor environment is 0.26. This is compared to between
% environment differences of max 0.51 (for low yield patches)

% much higher difference between medium and high yield patches in poor
% environment - SD difference = 0.66. 

%% test simulated data 

% load simulated data 
load ~/Dropbox/foraging/outputs/M2/SimulationResults.mat

simSD_LT = zeros([numel(SimData), 3, 2]); % each row is subject SD, each column patch type, 3d is block type
for iSim = 1:numel(SimData) % for each simulation (participant best fit)
    for iPatch = 1:3
        for iBlock = 1:2
            PatchIndex = SimData{iSim}{iBlock}.PatchOrder == iPatch;
                simSD_LT(iSim,iPatch,iBlock) = std(SimData{iSim}{iBlock}.LeavingTime(PatchIndex));
        end
    end
end

% between environments 

sim_change_in_SD_between_environments = mean(simSD_LT(:,:,2) - simSD_LT(:,:,1)); % poor - rich
% differences in real SDs are again much smaller than analytical SDs, and don't exactly align with real participant SDs 

% within environments 

sim_change_in_SD_within_rich = mean(simSD_LT(:,:,1)) - mean(simSD_LT(:,:,1))';
% greater differences in SD between patches WITHIN rich environment, than between
% environments - as expected from analytical results 
sim_change_in_SD_within_poor = mean(simSD_LT(:,:,2)) - mean(simSD_LT(:,:,2))';
% but lower differences in SD between patches within poor environment,
% compared to between environments - contradicting analytical results but
% matching with real participant data