#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS,stats,tidyverse,brms)

# set seed
set.seed(42) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Import the data ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.interval.combined = read.csv("../data/Exp2DataInterval.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Change the variable names and contrast code ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EfficacyLag = c("Efficacy_T0",
                "Efficacy_T1",
                "Efficacy_T2",
                "Efficacy_T3",
                "Efficacy_T4")


df.interval.combined = df.interval.combined %>% 
  mutate_at(EfficacyLag,funs(recode(.,"0" = "NoEfficacy","1" = "Efficacy"))) %>% 
  mutate_at(EfficacyLag,funs(factor(.)))

# reorder the levels
df.interval.combined$Efficacy_T0 = ordered(df.interval.combined$Efficacy_T0, levels = c("NoEfficacy", "Efficacy"))
df.interval.combined$Efficacy_T1 = ordered(df.interval.combined$Efficacy_T1, levels = c("NoEfficacy", "Efficacy"))
df.interval.combined$Efficacy_T2 = ordered(df.interval.combined$Efficacy_T2, levels = c("NoEfficacy", "Efficacy"))
df.interval.combined$Efficacy_T3 = ordered(df.interval.combined$Efficacy_T3, levels = c("NoEfficacy", "Efficacy"))
df.interval.combined$Efficacy_T4 = ordered(df.interval.combined$Efficacy_T4, levels = c("NoEfficacy", "Efficacy"))

# contrast code
contrasts(df.interval.combined$Efficacy_T0) = contr.sdif(2)
contrasts(df.interval.combined$Efficacy_T1) = contr.sdif(2)
contrasts(df.interval.combined$Efficacy_T2) = contr.sdif(2)
contrasts(df.interval.combined$Efficacy_T3) = contr.sdif(2)
contrasts(df.interval.combined$Efficacy_T4) = contr.sdif(2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Set the priors###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prior = c(
  prior(normal(0.5,0.2), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 0.2), class = b)) # a wide prior

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Fit the efficacy response model ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

efficacy = brm(EfficacyProbeResp ~  1 +
                 Efficacy_T0 +
                 Efficacy_T1 +
                 Efficacy_T2 +
                 Efficacy_T3 +
                 Efficacy_T4 +
                 Reward_T0 +
                 Reward_T1 +
                 Reward_T2 +
                 Reward_T3 +
                 Reward_T4 +
                 (1|SubID),
               data=df.interval.combined[!is.na(df.interval.combined$EfficacyProbeResp),],
               family=gaussian(),
               prior = prior,
               iter  = 20000,
               warmup = 19000,
               save_pars=save_pars(all=TRUE),
               control = list(adapt_delta = 0.99),
               sample_prior = TRUE,
               cores=4)
saveRDS(efficacy,file="../output/model.efficacy.previous.feedback.rds")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Fit the reward response model ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Predicting reward probe response by reward and efficacy feedback
reward = brm(RewardProbeResp ~  1 +
               Reward_T0 +
               Reward_T1 +
               Reward_T2 +
               Reward_T3 +
               Reward_T4 +
               Efficacy_T0 +
               Efficacy_T1 +
               Efficacy_T2 +
               Efficacy_T3 +
               Efficacy_T4 +
               (1|SubID),
             data=df.interval.combined[!is.na(df.interval.combined$RewardProbeResp),],
             family=gaussian(),
             prior = prior,
             iter  = 20000,
             warmup = 19000,
             save_pars=save_pars(all=TRUE),
             control = list(adapt_delta = 0.99),
             sample_prior = TRUE,
             cores=4)
saveRDS(reward,file="../output/model.reward.previous.feedback.rds")