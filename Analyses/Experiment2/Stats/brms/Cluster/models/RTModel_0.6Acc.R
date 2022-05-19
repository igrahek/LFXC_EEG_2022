library("brms")

# set seed
set.seed(42) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Import the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.trial.combined = read.csv("../data/Exp2Data_0.6Acc.csv")

# Change the name of the congruency variable
df.trial.combined$Congruency = df.trial.combined$type

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Set the priors based on Experiment 1 posteriors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior = c(
  prior(normal(655.36, 10.21), class = Intercept),
  prior(normal(71.29, 5.81), class = b, coef=Congruencyincongruent),
  prior(normal(-16.00, 9.47), class = b, coef=mbased_efficacy_prev),
  prior(normal(-1.56, 13.32), class = b, coef=mbased_reward_prev))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Fit the RT model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.congruency_plus_efficacy.RT = brm(AccRT ~ Congruency + mbased_efficacy_prev + mbased_reward_prev + scaledIntervalLength +  scaledIntervalCong+ (Congruency + mbased_efficacy_prev + mbased_reward_prev|SubID),
                                        data=df.trial.combined[!is.na(df.trial.combined$mbased_reward_prev) & !is.na(df.trial.combined$AccRT),],
                                        family=exgaussian(),
                                        prior = prior,
                                        cores = 4,
                                        iter  = 20000,
                                        warmup = 19000,
                                        save_pars=save_pars(all=TRUE),
                                        sample_prior = TRUE)
saveRDS(model.congruency_plus_efficacy.RT,file="../output/model.congruency_plus_efficacy.RT_0.7Acc.rds")
