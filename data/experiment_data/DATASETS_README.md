# README: Patch foraging datasets

## Files and descriptions: 

We have access to three datasets: Le Heron et al (2020), Contreras-Huerta et al (2022), and Kane et al (2019) (publicly available data from Figure 1 of their paper). 

There are 4 derived files associated with each dataset:

 ### trialbytrial.csv
A CSV file containing the leaving times for each subject, for each individual patch. The key variables used are:
- sub = subject
- patch = patch type [1 = low, 2 = medium, 3 = high yield]
- env = environment type [1 = rich, 2 = poor]
- leaveT = seconds stayed in patch (leaving time) or number of harvest in patch (rat data)
- [Note: The conterashuerta file has an additional column ‘ben’, which is ‘beneficiary’ [1 = collecting rewards for self, 2 = collecting rewards for another person]. Our paper only analyses data for the self condition (ben == 1).]

### subj_var.mat
A MATLAB struct variable with further information relevant to each subject required for model fitting:

- blockOrder - an array of the counterbalanced block order [1 = rich environment, 2 = poor environment]. Each row is a subject, each column is a block.

- blockSwitchIndex - a cell array (n = number of subjects). Each cell contains an index for when the block switches in the task for that subject. Each row is a new patch. 1 = block switch occurs. 

- experiencedAvgRR - an array containing the average reward rate achieved for each subject, across each type of environment. First column = rich environment, second column = poor environment. Each row is a subject. 

### CV_early_late.csv
A CSV file containing the calculated coefficient of variation for each subject in early [column 1] vs late patches in the task [column 2]. For further analysis in JASP.

### subject_LT.mat
A MATLAB structure containing the average leaving time and standard deviation of leaving time for each subject, for each patch and environment. Dimensions are [n_environment, n_patch, n_subject]. 



