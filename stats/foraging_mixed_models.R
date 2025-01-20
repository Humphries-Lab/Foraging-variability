#-- Mixed models for foraging variability paper - Figure 1 and 5 --# 
# Emma Scholey, 4 Nov 2024
# update 15 January 2025 - 

rm(list = ls(all = TRUE)) # clear environment

if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS, tidyverse, lme4, ggpubr, BayesFactor)

study = 'contrerashuerta' # v1, v3 or mri. 

#read data
data <- read.csv(paste('/Users/exs165/Dropbox/foraging/foraging_variability/data/experiment_data/',study,'/',study,'_trialbytrial.csv', sep = ""))

n_patch = 3

if (study == 'contrerashuerta'){
  # Filter out other beneficiary in Contreras-huerta dataset
  data <- data %>% filter(ben == 1)
  n_patch = 2
}

# MIXED MODELS ##############
data$mean_lt <- data$leaveT

##--------------------------- prepare/recode data for mixed models -----------------------
#
# # plots - check still looks okay
# # ggplot(data, aes(x = patch, y = mean_lt, group = env)) +
# #   stat_summary(fun.data = "mean_se", geom = "pointrange", size = 0.8)
#
# # check distribution
# ggplot(data, aes(mean_lt)) +
#   geom_histogram() + facet_wrap(~ patch + env)
#
# # set contrasts for patch and environment
data$patch <- factor(data$patch, ordered = T)
data$env <- factor(data$env, ordered = T)

contrasts(data$patch) <- contr.sdif(n_patch) # repeated contrasts ('sliding differences'). CHANGE to 2 for Contreras-Huerta
contrasts(data$env) <- contr.sdif(2) # same as coding block type as 0.5 or -0.5

# ##--------------------------- run main mixed model -----------------------
# ## Means - Figure 1

m_mean <- lmer(mean_lt ~ patch*env +
                    (patch*env|sub), data = data,
                 control = lmerControl(optimizer = c('bobyqa'), calc.derivs = TRUE, optCtrl=list(maxfun=20000)))
summary(m_mean)

#
if (n_patch == 3){
  tmp <- model.matrix(m_mean)
  patch2_1 <- model.matrix(m_mean)[,2]
  patch3_2 <- model.matrix(m_mean)[,3]
  env2_1 <- model.matrix(m_mean)[,4]
  patch2_1_env <- model.matrix(m_mean)[,5]
  patch3_2_env <- model.matrix(m_mean)[,6]
} else if (n_patch == 2){
  tmp <- model.matrix(m_mean)
  patch2_1 <- model.matrix(m_mean)[,2]
  env2_1 <- model.matrix(m_mean)[,3]
  patch2_1_env <- model.matrix(m_mean)[,4]
  
}
#

if (study == 'leheron'){
# # LE HERON
m_mean <- lmer(mean_lt ~ patch*env +
             (patch2_1 + patch3_2 + env2_1 |sub), data = data,
           control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_mean)


# get statistics for main effects 
m_mean <- lmerTest::lmer(mean_lt ~ patch * env + 
                           (patch2_1 + patch3_2 + env2_1|sub), data = data,
                         control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
anova(m_mean)
summary(m_mean)


} else if (study == 'contrerashuerta'){
# CONTRERAS-HUERTA
m_mean <- lmer(mean_lt ~ patch*env +
                 (patch*env|sub), data = data,
               control = lmerControl(optimizer = c('bobyqa'), calc.derivs = TRUE, optCtrl=list(maxfun=20000)))
summary(m_mean)

m_mean <- lmerTest::lmer(mean_lt ~ patch*env +
                            (patch*env|sub), data = data,
                          control = lmerControl(optimizer = c('bobyqa'), calc.derivs = TRUE, optCtrl=list(maxfun=20000)))
anova(m_mean)
summary(m_mean)


} else if (study == 'kane'){
# # KANE
m_mean <- lmer(mean_lt ~ patch*env +
             (patch2_1 + env2_1 + patch3_2_env + patch2_1_env||sub), data = data,
           control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_mean)

m_mean <- lmerTest::lmer(mean_lt ~ patch*env +
                       (patch2_1 + env2_1 + patch3_2_env + patch2_1_env||sub), data = data,
                     control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
anova(m_mean)
summary(m_mean)
}

# #------ model diagnostics
# 
# sum(abs(resid(m_mean, scaled = TRUE)) > 2) / length(resid(m_mean)) # should be max 5%
# sum(abs(resid(m_mean, scaled = TRUE)) > 2.5) / length(resid(m_mean)) # should be max 1%
# sum(abs(resid(m_mean, scaled = TRUE)) > 3) / length(resid(m_mean)) # should be 0%
# 
# # autocorrelation
# plot(acf(data$mean_lt))
# plot(acf(resid(m_mean)))
# 
# # homoscedasticity
# plot(m_mean,type=c('p','smooth'))
# 

# ##--------------------------- run SD leave model ----------------------- 
data <- data %>% group_by(sub, patch, env) %>% summarise(sd_lt = sd(leaveT))

ggplot(data, aes(sd_lt)) + 
  geom_histogram() + facet_wrap(~ patch + env)

# set contrasts for patch and environment
data$patch <- factor(data$patch, ordered = T)
data$env <- factor(data$env, ordered = T)

contrasts(data$patch) <- contr.sdif(n_patch) # repeated contrasts ('sliding differences')
contrasts(data$env) <- contr.sdif(2) # same as coding block type as 0.5 or -0.5

## Standard deviations - Figure 5
m_sd <- lmer(sd_lt ~ patch*env +
             (patch+env|sub), data = data,
           control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_sd)

if (n_patch == 3){
  tmp <- model.matrix(m_sd)
  patch2_1 <- model.matrix(m_sd)[,2]
  patch3_2 <- model.matrix(m_sd)[,3]
  env2_1 <- model.matrix(m_sd)[,4]
  patch2_1_env <- model.matrix(m_sd)[,5]
  patch3_2_env <- model.matrix(m_sd)[,6]
} else if (n_patch == 2){
  tmp <- model.matrix(m_sd)
  patch2_1 <- model.matrix(m_sd)[,2]
  env2_1 <- model.matrix(m_sd)[,3]
  patch2_1_env <- model.matrix(m_sd)[,4]
  
}

## LE HERON

if (study=='leheron'){
m_sd <- lmer(sd_lt ~ patch*env +
             (0 + patch2_1 + patch3_2 + env + patch2_1_env + patch3_2_env||sub), data = data,
           control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_sd)

m_sd <- lmerTest::lmer(sd_lt ~ patch*env +
                       (0 + patch2_1 + patch3_2 + env + patch2_1_env + patch3_2_env||sub), data = data,
                     control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
anova(m_sd)
summary(m_sd)

data$sub <- factor(data$sub)

# BF model with 
bf_full = lmBF(sd_lt ~ patch*env + sub, data = data, whichRandom = c('sub'))

# BF model without 
bf_no_patch = lmBF(sd_lt ~ env + patch:env + sub, data = data, whichRandom = c('sub'))
bf_no_env = lmBF(sd_lt ~ patch + patch:env + sub, data = data, whichRandom = c('sub'))

bf_patch = 1/(bf_full/bf_no_patch)
bf_env = 1/(bf_full/bf_no_env)

bf_patch
bf_env
} else if (study =='contrerashuerta'){
## CONTRERAS-HUERTA

m_sd <- lmer(sd_lt ~ patch*env +
               (patch2_1 + env2_1||sub), data = data,
             control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_sd)

m_sd <- lmerTest::lmer(sd_lt ~ patch*env +
                         (patch2_1 + env2_1||sub), data = data,
                       control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
anova(m_sd)
summary(m_sd)

data$sub <- factor(data$sub)

# BF model with 
bf_full = lmBF(sd_lt ~ patch*env + sub, data = data, whichRandom = c('sub', 'patch', 'env'))

# BF model without 
bf_no_patch = lmBF(sd_lt ~ env + patch:env + sub, data = data, whichRandom = c('sub', 'patch', 'env'))
bf_no_env = lmBF(sd_lt ~ patch + patch:env + sub, data = data, whichRandom = c('sub', 'patch', 'env'))

bf_patch = 1/(bf_full/bf_no_patch)
bf_env = 1/(bf_full/bf_no_env)

bf_patch
bf_env

} else if (study=='kane'){
## KANE
m_sd <- lmer(sd_lt ~ patch*env +
               (env2_1 + patch2_1_env + patch3_2_env||sub), data = data,
             control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
summary(m_sd)

m_sd <- lmerTest::lmer(sd_lt ~ patch*env +
                         (env2_1 + patch2_1_env + patch3_2_env||sub), data = data,
                       control = lmerControl(optimizer = c('bobyqa'), optCtrl=list(maxfun=20000)))
anova(m_sd)
summary(m_sd)

data$sub <- factor(data$sub)

# BF model with patch
bf_full = lmBF(sd_lt ~ patch*env + sub, data = data, whichRandom = c('sub'))

# BF model without 
bf_no_patch = lmBF(sd_lt ~ env + patch:env + sub, data = data, whichRandom = c('sub'))
bf_no_env = lmBF(sd_lt ~ patch + patch:env + sub, data = data, whichRandom = c('sub'))

bf_patch = 1/(bf_full/bf_no_patch)
bf_env = 1/(bf_full/bf_no_env)

bf_patch
bf_env

}