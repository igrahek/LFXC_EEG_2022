#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(MASS,stats,tidyverse,brms)

# set seed
set.seed(42) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Import the data ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
learning_rates_efficacy = read.csv("../data/Exp2LRsEfficacy.csv")
learning_rates_reward = read.csv("../data/Exp2LRsReward.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Efficacy ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# wide to long format
lr_ttest = subset(learning_rates_efficacy, select=c("positive_learning_rate","negative_learning_rate","subject"))
lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)

# set the variable types
lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)

# Set the priors
prior = c(
  prior(normal(0.5,0.5), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 0.5), class = b)) # a wide prior

# All subjects
m = brm(
  Value ~ Learning_rate_type, sigma ~ Learning_rate_type,
  prior=prior,
  family=student,
  data = lr_ttest,
  iter  = 20000,
  warmup = 19000,
  save_pars=save_pars(all=TRUE),
  control = list(adapt_delta = 0.99),
  sample_prior = TRUE,
  cores=4)
saveRDS(m,file="../output/model.ttest.learning_rates_efficacy.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Reward ###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# wide to long format
lr_ttest = subset(learning_rates_reward, select=c("positive_learning_rate","negative_learning_rate","subject"))
lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)

# set the variable types
lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)

# Set the priors
prior = c(
  prior(normal(0.5,0.5), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 0.5), class = b)) # a wide prior

m = brm(
  Value ~ Learning_rate_type, sigma ~ Learning_rate_type,
  prior=prior,
  family=student,
  data = lr_ttest,
  iter  = 20000,
  warmup = 19000,
  cores=4,
  save_pars=save_pars(all=TRUE),
  sample_prior = TRUE)
saveRDS(m,file="../output/model.ttest.learning_rates_reward.rds")
