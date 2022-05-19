library("brms")

# set seed
set.seed(42) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Import the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.trial.combined = read.csv("../data/Exp2Data_0.6Acc.csv")

df.trial.combined$hit = as.numeric(df.trial.combined$hit)

# Change the name of the congruency variable
df.trial.combined$Congruency = df.trial.combined$type

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Set the priors based on Experiment 1 posteriors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior = c(
  prior(normal(1.56, 0.12), class = Intercept),
  prior(normal(-0.59, 0.09), class = b, coef=Congruencyincongruent),
  prior(normal(0.11, 0.16), class = b, coef=mbased_efficacy_prev),
  prior(normal(-0.36, 0.19), class = b, coef=mbased_reward_prev))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Fit the accuracy model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.congruency_plus_efficacy.Acc = brm(hit ~ Congruency + mbased_efficacy_prev + mbased_reward_prev + scaledIntervalLength +  scaledIntervalCong+  (Congruency + mbased_efficacy_prev + mbased_reward_prev|SubID),
                                         data=df.trial.combined[!is.na(df.trial.combined$response) & !is.na(df.trial.combined$mbased_reward_prev),],
                                         family="bernoulli",
                                         prior = prior,
                                         iter  = 20000,
                                         warmup = 19000,
                                         save_pars=save_pars(all=TRUE),
                                         cores = 4,
                                         sample_prior = TRUE)
saveRDS(model.congruency_plus_efficacy.Acc,file="../output/model.congruency_plus_efficacy.Acc_0.7Acc.rds")
