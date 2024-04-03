#-- Standard deviation equivalence testing for stochastic foraging paper --# 
# To test that SDs are equivalent using the TOST method 

# RESOURCES 
# https://github.com/Lakens/TOSTER?tab=readme-ov-file
# https://aaroncaldwell.us/TOSTERpkg/index.html
# https://cloud.r-project.org/web/packages/TOSTER/TOSTER.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5502906/#bibr18-1948550617697177 - paper discussing TOST 
# https://journals.sagepub.com/doi/10.1177/2515245918770963 - principles for setting equivalence bounds based on SESOI
# Emma Scholey
# 20 February 2024


library(TOSTER)
library(tidyverse)
library(BayesFactor)
library(rstatix)
library(ggpubr)
library(emmeans)
library(kableExtra)

options(knitr.table.format = "latex") 
#------------------------- early vs late SD ----------------------------------# 

#################### Le Heron - early vs late SD 

d_wide <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig2_leheron_SD_early_late_data.csv', header = FALSE)
colnames(d_wide) <- c('early', 'late')
d <- d_wide %>% pivot_longer(cols = c(1,2), names_to = 'phase')

# calculate difference scores 
diffSD <- with(d, value[phase == 'early'] - value[phase == 'late'])
# check normality of differences
shapiro_test(diffSD) # not significantly different, assume normality 

# based on SD of differences, compute the minimal mean raw diff to find a statistically significant effect (based on minimum effect size)
# otherwise known as SESOI - smallest effect size of interest
# minimum effect size is G*POWER calculations, to find Cohen's d for this sample size, alpha = .05 and power = 0.80. 
sdDiff <- sd(diffSD)
minimal_meanDiff = sdDiff*0.46 # since Cohen's d = meanDiff/stdDiff

# run equivalence test 
nulltest_raw <- t_TOST(formula = value ~ phase, data = d, paired = TRUE, eqb = minimal_meanDiff)
describe(nulltest_raw)

# compare to Bayes Factors

## Traditional two-tailed t test
t.test(diffSD)
bf = ttestBF(x = d_wide$early,d_wide$late, paired = TRUE)
bf

################# Contreras Huerta - early vs late SD 

d_wide <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig2_contrerashuerta_SD_early_late_data.csv', header = FALSE)
colnames(d_wide) <- c('early', 'late')
d <- d_wide %>% pivot_longer(cols = c(1,2), names_to = 'phase')

# calculate difference scores 
diffSD <- with(d, value[phase == 'early'] - value[phase == 'late'])
# check normality of differences
shapiro_test(diffSD) # not significantly different, assume normality 

# based on SD of differences, compute the minimal mean raw diff to find a statistically significant effect (based on minimum effect size)
# otherwise known as SESOI - smallest effect size of interest
# minimum effect size is G*POWER calculations, to find Cohen's d for this sample size, alpha = .05 and power = 0.80. 
sdDiff <- sd(diffSD)
minimal_meanDiff = sdDiff*0.54 

# run equivalence test 
nulltest_raw <- t_TOST(formula = value ~ phase, data = d, paired = TRUE, eqb = minimal_meanDiff)
describe(nulltest_raw)

# run Bayes Factor
t.test(diffSD)
bf = ttestBF(x = d_wide$early,d_wide$late, paired = TRUE)
bf

###################### Kane- early vs late SD 

d_wide <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig2_kane_SD_early_late_data.csv', header = FALSE)
colnames(d_wide) <- c('early', 'late')
d <- d_wide %>% pivot_longer(cols = c(1,2), names_to = 'phase')

# calculate difference scores 
diffSD <- with(d, value[phase == 'early'] - value[phase == 'late'])
# check normality of differences
shapiro_test(diffSD) # not significantly different, assume normality 

# based on SD of differences, compute the minimal mean raw diff to find a statistically significant effect (based on minimum effect size)
# otherwise known as SESOI - smallest effect size of interest
# minimum effect size is G*POWER calculations, to find Cohen's d for this sample size, alpha = .05 and power = 0.80. 
sdDiff <- sd(diffSD)
minimal_meanDiff = sdDiff*1.16 

# run equivalence test 
nulltest_raw <- t_TOST(formula = value ~ phase, data = d, paired = TRUE, eqb = minimal_meanDiff)
describe(nulltest_raw)

# run Bayes Factor
t.test(diffSD)
bf = ttestBF(x = d_wide$early,d_wide$late, paired = TRUE)
bf


#------------------------- SD between env/patch ----------------------------------# 

################## --------------------- Le Heron 

##---------------- Empirical data
d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_leheron_subject_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','rich_high','poor_low','poor_mid','poor_high')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:6), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)

# check assumptions of anova
# normality 
d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# not normal for most cells - long right tail (since working with timings) - use logs 
d <- d %>% mutate(log_value = log(value))
d %>% group_by(env,patch) %>% shapiro_test(log_value)
ggqqplot(d, "log_value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(value,type = "mean_sd")

# run emmeans equivalence test 
afmod <- afex::aov_car(log_value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = 'bonf')
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor

##------------------ Simulated data

d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_leheron_sim_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','rich_high','poor_low','poor_mid','poor_high')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:6), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)

# check assumptions of anova
# normality 
d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# not normal for most cells - long right tail (since working with timings) - use logs 
d <- d %>% mutate(log_value = log(value))
d %>% group_by(env,patch) %>% shapiro_test(log_value)
ggqqplot(d, "log_value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(value,type = "mean_sd")

# run emmeans equivalence test 
afmod <- afex::aov_car(log_value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = "bonf")
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor


################# ---------------------------------Contreras-Huerta

##---------------- Empirical data
d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_contrerashuerta_subject_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','poor_low','poor_mid')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:4), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)


# check assumptions of anova
# normality 
d <- d %>% mutate(log_value = log(value))

d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# not normal - long right tail (since working with timings) - use logs 
d %>% group_by(env,patch) %>% shapiro_test(log_value)
ggqqplot(d, "log_value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(log_value,type = "mean_sd")

# run ANOVA 
afmod <- afex::aov_car(log_value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = 'bonf')
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor

##---------------------- Simulated data
d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_contrerashuerta_sim_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','poor_low','poor_mid')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:4), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)


# check assumptions of anova
# normality 
d <- d %>% mutate(log_value = log(value))

d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# not normal - long right tail (since working with timings) - use logs 
d %>% group_by(env,patch) %>% shapiro_test(log_value)
ggqqplot(d, "log_value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(log_value,type = "mean_sd")

# run ANOVA 
afmod <- afex::aov_car(log_value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = 'bonf')
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor


#################---------------------------------- Kane 

##---------------------- Empirical data
d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_kane_subject_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','rich_high','poor_low','poor_mid','poor_high')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:6), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)


# check assumptions of anova
# normality 
d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# normal - use value (not log)

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(value,type = "mean_sd")


# run emmeans equivalence test 
afmod <- afex::aov_car(value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = 'bonf')
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor

##------------------------ Simulated data
d <- read.csv('/Users/exs165/Dropbox/foraging/stochastic-mvt-project/data/experiment_data/figures_stats_data/fig4_kane_sim_SD_data.csv', header = FALSE)
colnames(d) <- c('rich_low','rich_mid','rich_high','poor_low','poor_mid','poor_high')

d$subj <- 1:nrow(d) # add subject number to include as random effect
d <- d %>% pivot_longer(cols = c(1:6), names_to = c('env', 'patch'), names_sep = '_')
d$env <- factor(d$env)
d$patch <- factor(d$patch, levels = c('low', 'mid','high'))
d$subj <- factor(d$subj)


# check assumptions of anova
# normality 
d %>% group_by(env,patch) %>% shapiro_test(value)
ggqqplot(d, "value", ggtheme = theme_bw()) + facet_grid(patch ~ env, labeller = "label_both")
# normal - use value (not log)

# report summary statistics 
d %>% group_by(env,patch) %>% get_summary_stats(value,type = "mean_sd")


# run emmeans equivalence test 
afmod <- afex::aov_car(value ~ env*patch + Error(subj/(env*patch)), data = d, anova_table = list(es = "pes"))
afmod
afex::nice(afmod) %>% kbl(format = 'latex') %>% kable_classic()

# marginal means
em_patch <- emmeans::emmeans(afmod, specs = ~ patch)
em_env <- emmeans::emmeans(afmod, specs = ~ env)
em_interaction <- emmeans::emmeans(afmod, specs = ~ patch*env)

# contrasts - difference between marginal means 
pairs(em_patch, adjust = 'bonf')
pairs(em_env, adjust = 'bonf')
pairs(em_interaction, adjust = 'bonf')
plot(em_patch)

# run Bayes Factor ANOVA 
bf = anovaBF(formula = log_value ~ env*patch + subj, data = d, whichRandom = "subj")
bf
bf[4] / bf[3] # calculate interaction Bayes Factor