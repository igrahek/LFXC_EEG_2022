#### About the code ####

# Experiment: LFXC_EEG 
# Code written by: Ivan Grahek (2018-2021)
# Description: Code for the analysis of behavioral data for the LFXC_EEG project.  

#### Clear the environment and import the data ####---------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan,plyr,Rmisc,yarrr,BayesFactor,reshape2,brms, broom, tidyverse, knitr, here, zoo, magrittr, pracma,xtable, Hmisc, ppcor, lme4,MuMIn,MASS, sjPlot, jtools, lmerTest, sjstats,coefplot,R.matlab,RColorBrewer,ggeffects,glmmTMB, interactions)

# stop the build tools pop up
options(buildtools.check = function(action) TRUE )

# set seed
set.seed(42) 
# set directory
setwd(here())

# import data
data.raw = read.csv(file = here("Data/Experiment1/Behavior/preprocessed","data104.csv"),header=TRUE,na.strings="NaN")

#rename variables and save as factors
data.raw <- data.raw %>% 
  mutate_at(c("Congruency"), funs(recode(., "0" = "Incongruent","1" = "Congruent","2" = "Neutral")))

RewardLag = c("IsRewarded_T1",
              "IsRewarded_T2",
              "IsRewarded_T3",
              "IsRewarded_T4",
              "IsRewarded_T5",
              "IsRewarded_T6",
              "IsRewarded_T7",
              "IsRewarded_T8",
              "IsRewarded_T9",
              "IsRewarded_T10")

EfficacyLag = c("Efficacy_T1",
                "Efficacy_T2",
                "Efficacy_T3",
                "Efficacy_T4",
                "Efficacy_T5",
                "Efficacy_T6",
                "Efficacy_T7",
                "Efficacy_T8",
                "Efficacy_T9",
                "Efficacy_T10")



# Recode and turn into a factor
data.raw = data.raw %>% 
  mutate_at(RewardLag,funs(recode(.,"0" = "NoReward","1" = "Reward"))) %>% 
  mutate_at(RewardLag,funs(factor(.)))

data.raw = data.raw %>% 
  mutate_at("IsRewarded",funs(recode(.,"0" = "NoReward","1" = "Reward"))) %>% 
  mutate_at("IsRewarded",funs(factor(.)))

data.raw = data.raw %>% 
  mutate_at(EfficacyLag,funs(recode(.,"0" = "NoEfficacy","1" = "Efficacy"))) %>% 
  mutate_at(EfficacyLag,funs(factor(.)))

data.raw = data.raw %>% 
  mutate_at("EffLvl",funs(recode(.,"0" = "NoEfficacy","1" = "Efficacy"))) %>% 
  mutate_at("EffLvl",funs(factor(.)))


# Create contrasts for the n-back regression plots
data.raw$IsRewarded = ordered(data.raw$IsRewarded, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T1 = ordered(data.raw$IsRewarded_T1, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T2 = ordered(data.raw$IsRewarded_T2, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T3 = ordered(data.raw$IsRewarded_T3, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T4 = ordered(data.raw$IsRewarded_T4, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T5 = ordered(data.raw$IsRewarded_T5, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T6 = ordered(data.raw$IsRewarded_T6, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T7 = ordered(data.raw$IsRewarded_T7, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T8 = ordered(data.raw$IsRewarded_T8, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T9 = ordered(data.raw$IsRewarded_T9, levels = c("NoReward", "Reward"))
data.raw$IsRewarded_T10 = ordered(data.raw$IsRewarded_T10, levels = c("NoReward", "Reward"))

data.raw$EffLvl = ordered(data.raw$EffLvl, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T1 = ordered(data.raw$Efficacy_T1, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T2 = ordered(data.raw$Efficacy_T2, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T3 = ordered(data.raw$Efficacy_T3, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T4 = ordered(data.raw$Efficacy_T4, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T5 = ordered(data.raw$Efficacy_T5, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T6 = ordered(data.raw$Efficacy_T6, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T7 = ordered(data.raw$Efficacy_T7, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T8 = ordered(data.raw$Efficacy_T8, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T9 = ordered(data.raw$Efficacy_T9, levels = c("NoEfficacy", "Efficacy"))
data.raw$Efficacy_T10 = ordered(data.raw$Efficacy_T10, levels = c("NoEfficacy", "Efficacy"))


data.raw = ddply(data.raw,.(SubID),plyr::mutate,
                 Trial = 1:(length(data.raw$RT)/length(unique(data.raw$SubID))))

# create a new variable to account for different drifts that subjects have
data.raw$Drift = ifelse(data.raw$SubID %% 2 == 1, "Odd","Even")
data.raw$Drift = as.factor(data.raw$Drift)

# NA the RTs that have been identified to be with wrong latencies
data.raw$RT[4607] = NA
data.raw$RT[8090] = NA
data.raw$RT[8258] = NA

# # Create the variable with three levels to analyze efficacy and reward feedback taking into account the feedback order
data.raw$Rew_eff_fb[data.raw$feedbackOrder == 1]="RewardUnknown"
data.raw$Rew_eff_fb[data.raw$feedbackOrder == 2 & data.raw$IsRewarded == "Reward"]="Reward"
data.raw$Rew_eff_fb[data.raw$feedbackOrder == 2 & data.raw$IsRewarded == "NoReward"]="NoReward"

data.raw$Eff_rew_fb[data.raw$feedbackOrder == 2]="EfficacyUnknown"
data.raw$Eff_rew_fb[data.raw$feedbackOrder == 1 & data.raw$EffLvl == "Efficacy"]="Efficacy"
data.raw$Eff_rew_fb[data.raw$feedbackOrder == 1 & data.raw$EffLvl == "NoEfficacy"]="NoEfficacy"

## Add the error type variable

## compute 3-level variable: correct, distractor error, random error

data.raw$Distractorerror <- rep(0, length(data.raw$SubID))
data.raw$Distractorerror[data.raw$Acc==0 & data.raw$Resp==data.raw$DistractorResp] <- 1

data.raw$Randomerror <- rep(0, length(data.raw$SubID))
data.raw$Randomerror[data.raw$Acc==0 & !data.raw$Resp==data.raw$DistractorResp] <- 1

data.raw$Errortype <- rep("correct", length(data.raw$SubID))
data.raw$Errortype[data.raw$Distractorerror==1] <- "distractor error"
data.raw$Errortype[data.raw$Randomerror==1] <- "random error"
data.raw$Errortype[data.raw$isMiss==1] <- "miss"

dim(subset(data.raw,data.raw$ErrorType=="Automatic"))
# ggplot(data= data.raw, aes(x = Errortype)) + geom_bar()+ theme_classic()



#### Save data for the RL model fitting ####

# Edit the dataset
data_RLfit = data.raw %>%
  dplyr::select(SubID, 
                EffLvl,
                IsRewarded,
                EfficacyProbeResp,
                RewRateProbeResp) 

# Change efficacy to 0-1
data_RLfit$EffLvl = ifelse(data_RLfit$EffLvl=="Efficacy",1,0)

# Change reward to 0-1
data_RLfit$IsRewarded = ifelse(data_RLfit$IsRewarded=="Reward",1,0)

data_RLfit = as.matrix(data_RLfit)


# Save out the data
write.csv(x = data_RLfit, file = paste0('Analyses/Experiment1/RL_fitting/Stan/Cluster/data/','RL_fit_data_','Experiment1','.csv'), row.names = FALSE)


#### Import the learning rates for efficacy ###### 
learning_rates = read.csv(file = here("Analyses/Experiment1/RL_fitting/Hierarchical/4_Intercept_Pos_and_neg_learning_rate/results","learning_rates_efficacy.csv"),header=F,na.strings="NaN")

# add the subject names
learning_rates$SubID = unique(data.raw$SubID)

# rename the variable names
colnames(learning_rates)[1] <- "positive_learning_rate"
colnames(learning_rates)[2] <- "negative_learning_rate"
colnames(learning_rates)[3] <- "initial_bias"
colnames(learning_rates)[4] <- "subject"
colnames(learning_rates)[5] <- "BIC"

learning_rates_efficacy = learning_rates

# add a new variable in the main data set in which efficacy is coded as 1 or 0
data.raw$EffLvl_forLR = ifelse(data.raw$EffLvl=="Efficacy",1,0)

#### Calculate the model-based efficacy estimate ######---------------------------------------------------------------------------------------------------------------------------------

# for each subject
for (s in 1:length(unique(data.raw$SubID))) {
  # for each trial
  for (t in 1:length(unique(data.raw$Trial))) {
    # if this is the first trial use the initial estimate
    if (t == 1) {
      v = learning_rates$initial_bias[s]
    } else {
      v = v
    }
    
    # save the BIC value
    data.raw$BIC_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                            data.raw$Trial == t] = learning_rates$BIC[s]
    
    # calculate the signed prediction error for this subject for this trial (difference between the actual and the expected efficacy)
    data.raw$signed_PE_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                                  data.raw$Trial == t] = data.raw$EffLvl_forLR[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                 data.raw$Trial == t] - v
    
    # calculate the unsigned prediction error for this subject for this trial (absolute difference between the actual and the expected efficacy)
    data.raw$unsigned_PE_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                                    data.raw$Trial == t] = abs(data.raw$EffLvl_forLR[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                       data.raw$Trial == t] - v)
    
    # save the difference between the positive and the negative LR
    data.raw$LR_efficacy_Pos_min_Neg[data.raw$SubID == unique(data.raw$SubID)[s] &
                                       data.raw$Trial == t] = learning_rates$positive_learning_rate[s] - learning_rates$negative_learning_rate[s]
    
    # if the signed pe is positive
    if (data.raw$signed_PE_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                                    data.raw$Trial == t] > 0) {
      
      
      
      # calculate the size of the update (the learning rate times the prediction error)
      data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                        data.raw$Trial == t] = learning_rates$positive_learning_rate[s] * data.raw$signed_PE_efficacy[data.raw$SubID ==
                                                                                                                                        unique(data.raw$SubID)[s] & data.raw$Trial == t]
      # save the absolute value of the update
      data.raw$mbased_efficacy_absolute_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                 data.raw$Trial == t] =  abs(data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                                               data.raw$Trial == t])
      
      # calculate the effiacy estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      data.raw$mbased_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                                 data.raw$Trial == t] = v + data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # update the value
      v = v + data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # save the learning rate value
      data.raw$LR_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                             data.raw$Trial == t] = learning_rates$positive_learning_rate[s]
      
      # save the positive learning rate value
      data.raw$LR_efficacy_positive[data.raw$SubID == unique(data.raw$SubID)[s] &
                                      data.raw$Trial == t] = learning_rates$positive_learning_rate[s]
      
      # if the rpe is negative
    } else {
      
      # calculate the size of the update (the learning rate times the prediction error)
      data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                        data.raw$Trial == t] = learning_rates$negative_learning_rate[s] * data.raw$signed_PE_efficacy[data.raw$SubID ==
                                                                                                                                        unique(data.raw$SubID)[s] & data.raw$Trial == t]
      # save the absolute value of the update
      data.raw$mbased_efficacy_absolute_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                 data.raw$Trial == t] =  abs(data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                                               data.raw$Trial == t])
      
      # calculate the effiacy estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      data.raw$mbased_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                                 data.raw$Trial == t] = v + data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # update the value
      v = v + data.raw$mbased_efficacy_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # save the learning rate value
      data.raw$LR_efficacy[data.raw$SubID == unique(data.raw$SubID)[s] &
                             data.raw$Trial == t] = learning_rates$negative_learning_rate[s]
      
      # save the positive learning rate value
      data.raw$LR_efficacy_negative[data.raw$SubID == unique(data.raw$SubID)[s] &
                                      data.raw$Trial == t] = learning_rates$negative_learning_rate[s]
      
    }
    
    
  }
}


#### Import the learning rates for reward rate######--------------------------------------------------------------------------------------------------------------------------------- 
learning_rates = read.csv(file = here("Analyses/Experiment1/RL_fitting/Hierarchical/4_Intercept_Pos_and_neg_learning_rate/results","learning_rates_reward.csv"),header=F,na.strings="NaN")

# add the subject names
learning_rates$SubID = unique(data.raw$SubID)

# rename the variable names
colnames(learning_rates)[1] <- "positive_learning_rate"
colnames(learning_rates)[2] <- "negative_learning_rate"
colnames(learning_rates)[3] <- "initial_bias"
colnames(learning_rates)[4] <- "subject"
colnames(learning_rates)[5] <- "BIC"

learning_rates_reward = learning_rates


# add a new variable in the main dataset in which efficacy is coded as 1 or 0
data.raw$IsRewarded_forLR = ifelse(data.raw$IsRewarded=="Reward",1,0)

#### Calculate the model-based reward rate estimate ######---------------------------------------------------------------------------------------------------------------------------------

# for each subject
for (s in 1:length(unique(data.raw$SubID))) {
  # for each trial
  for (t in 1:length(unique(data.raw$Trial))) {
    # if this is the first trial use the initial estimate
    if (t == 1) {
      v = learning_rates$initial_bias[s]
    } else {
      v = v
    }
    
    # save the BIC value
    data.raw$BIC_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                          data.raw$Trial == t] = learning_rates$BIC[s]
    
    # calculate the prediction error for this subject for this trial (difference between the actual and the expected efficacy)
    data.raw$signed_PE_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                                data.raw$Trial == t] = data.raw$IsRewarded_forLR[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                   data.raw$Trial == t] - v
    
    # calculate the unsigned prediction error for this subject for this trial (absolute difference between the actual and the expected efficacy)
    data.raw$unsigned_PE_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                                  data.raw$Trial == t] = abs(data.raw$IsRewarded_forLR[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                         data.raw$Trial == t] - v)
    
    # if the rpe is positive
    if (data.raw$signed_PE_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                                  data.raw$Trial == t] > 0) {
      
      # calculate the size of the update (the learning rate times the prediction error)
      data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                      data.raw$Trial == t] = learning_rates$positive_learning_rate[s] * data.raw$signed_PE_reward[data.raw$SubID ==
                                                                                                                                    unique(data.raw$SubID)[s] & data.raw$Trial == t]
      # save the absolute value of the update
      data.raw$mbased_reward_absolute_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                               data.raw$Trial == t] =  abs(data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                                           data.raw$Trial == t])
      
      # calculate the reward estimate (expected reward plus the update - the learning rate-weighted prediction error)
      data.raw$mbased_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                               data.raw$Trial == t] = v + data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # update the value
      v = v + data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # save the learning rate value
      data.raw$LR_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                           data.raw$Trial == t] = learning_rates$positive_learning_rate[s]
      
      # save the positive learning rate value
      data.raw$LR_reward_positive[data.raw$SubID == unique(data.raw$SubID)[s] &
                                    data.raw$Trial == t] = learning_rates$positive_learning_rate[s]
      
      # if the rpe is negative
    } else {
      
      # calculate the size of the update (the learning rate times the prediction error)
      data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                      data.raw$Trial == t] = learning_rates$negative_learning_rate[s] * data.raw$signed_PE_reward[data.raw$SubID ==
                                                                                                                                    unique(data.raw$SubID)[s] & data.raw$Trial == t]
      # save the absolute value of the update
      data.raw$mbased_reward_absolute_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                               data.raw$Trial == t] =  abs(data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] &
                                                                                                           data.raw$Trial == t])
      
      # calculate the reward estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      data.raw$mbased_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                               data.raw$Trial == t] = v + data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # update the value
      v = v + data.raw$mbased_reward_update[data.raw$SubID == unique(data.raw$SubID)[s] & data.raw$Trial == t]
      
      # save the learning rate value
      data.raw$LR_reward[data.raw$SubID == unique(data.raw$SubID)[s] &
                           data.raw$Trial == t] = learning_rates$negative_learning_rate[s]
      
      # save the positive learning rate value
      data.raw$LR_reward_negative[data.raw$SubID == unique(data.raw$SubID)[s] &
                                    data.raw$Trial == t] = learning_rates$negative_learning_rate[s]
      
    }
    
    
  }
}


# Create a previous value variable (for predicting everything before the new feedback participants rely on the value from the previous trial)
data.raw = ddply(data.raw,.(SubID),transform,
                 mbased_reward_prev = append(mbased_reward,NA,after=0)[-(length(mbased_reward)+1)],
                 mbased_efficacy_prev = append(mbased_efficacy,NA,after=0)[-(length(mbased_efficacy)+1)])

# create the variable for analyzing the relevant estimate at the feedback level (estimate is updated or not based on the feedback sequence)
data.raw$mbased_reward_feedback_locked = ifelse(data.raw$Rew_eff_fb == "RewardUnknown", data.raw$mbased_reward_prev,data.raw$mbased_reward)
data.raw$mbased_efficacy_feedback_locked = ifelse(data.raw$Eff_rew_fb == "EfficacyUnknown", data.raw$mbased_efficacy_prev,data.raw$mbased_efficacy)


#### Import the EEG data####---------------------------------------------------------------------------------------------------------------------------------

#### Late CNV (1000 - 1500ms after fix cross onset)
CNV10001500 = readMat(here("Data/Experiment1/EEG/preprocessed","CNV10001500.mat"))

CNV10001500 = as.data.frame(CNV10001500) # convert to data frame

colnames(CNV10001500)[c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                         16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                         46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                         61, 62, 63, 64, 65)]=
  
  c('Fp1', 'Fpz', 'Fp2', 'AF3', 'AFz', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8',
    'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T7', 'C5', 'C3', 'C1',
    'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10',
    'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO3', 'POz', 'PO4', 'O1', 'Oz', 'O2',
    'LO1', 'IO1', 'IO2', 'LO2', 'Cz')
# 
data.raw$CNV10001500 = apply(cbind (CNV10001500$Fz,CNV10001500$F1,CNV10001500$F2,
                                    CNV10001500$FCz,CNV10001500$FC1,CNV10001500$FC2,
                                    CNV10001500$Cz,CNV10001500$C1,CNV10001500$C2), 1, mean) # add the average activation in ROI to the behaviroal dataset



#### P3b to the reward feedback (350 - 500ms after feedback onset) 
P3bR350500 = readMat(here("Data/Experiment1/EEG/preprocessed","P3bR350500.mat"))


P3bR350500 = as.data.frame(P3bR350500) # convert to data frame

colnames(P3bR350500)[c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                        61, 62, 63, 64, 65)]=
  
  c('Fp1', 'Fpz', 'Fp2', 'AF3', 'AFz', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8',
    'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T7', 'C5', 'C3', 'C1',
    'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10',
    'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO3', 'POz', 'PO4', 'O1', 'Oz', 'O2',
    'LO1', 'IO1', 'IO2', 'LO2', 'Cz')


data.raw$P3bR350500 = apply(cbind (P3bR350500$Pz,P3bR350500$P1,P3bR350500$P2,
                                   P3bR350500$POz,P3bR350500$PO3,P3bR350500$PO4,
                                   P3bR350500$CPz,P3bR350500$CP1,P3bR350500$CP2), 1, mean) # add the average activation in ROI to the behaviroal dataset



#### P3b to the efficacy feedback (350 - 500ms after feedback onset) 
P3bE350500 = readMat(here("Data/Experiment1/EEG/preprocessed","P3bE350500.mat"))

P3bE350500 = as.data.frame(P3bE350500) # convert to data frame

colnames(P3bE350500)[c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                        61, 62, 63, 64, 65)]=
  
  c('Fp1', 'Fpz', 'Fp2', 'AF3', 'AFz', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8',
    'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T7', 'C5', 'C3', 'C1',
    'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10',
    'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO3', 'POz', 'PO4', 'O1', 'Oz', 'O2',
    'LO1', 'IO1', 'IO2', 'LO2', 'Cz')


data.raw$P3bE350500 = apply(cbind (P3bE350500$Pz,P3bE350500$P1,P3bE350500$P2,
                                   P3bE350500$POz,P3bE350500$PO3,P3bE350500$PO4,
                                   P3bE350500$CPz,P3bE350500$CP1,P3bE350500$CP2), 1, mean) # add the average activation in ROI to the behaviroal dataset

#### P3b to the efficacy feedback for the control analysis (200 - 300ms after feedback onset) 
P3bEcontrol200300 = readMat(here("Data/Experiment1/EEG/preprocessed","P3bEcontrol200300.mat"))



P3bEcontrol200300 = as.data.frame(P3bEcontrol200300) # convert to data frame

colnames(P3bEcontrol200300)[c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                        46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                        61, 62, 63, 64, 65)]=
  
  c('Fp1', 'Fpz', 'Fp2', 'AF3', 'AFz', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8',
    'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T7', 'C5', 'C3', 'C1',
    'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10',
    'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO3', 'POz', 'PO4', 'O1', 'Oz', 'O2',
    'LO1', 'IO1', 'IO2', 'LO2', 'Cz')


data.raw$P3bEcontrol200300 = apply(cbind (P3bEcontrol200300$Pz,P3bEcontrol200300$P1,P3bEcontrol200300$P2,
                                          P3bEcontrol200300$POz,P3bEcontrol200300$PO3,P3bEcontrol200300$PO4,
                                   P3bEcontrol200300$CPz,P3bEcontrol200300$CP1,P3bEcontrol200300$CP2), 1, mean) # add the average activation in ROI to the behavioral dataset



#### Clean the data ####---------------------------------------------------------------------------------------------------------------------------------

data = data.raw

# delete the RTs below 200ms
data$RT = ifelse(data$RT<200,NA,data$RT)

# Center the continuous variables 
data$Trial = scale(data$Trial/288, scale= FALSE, center = TRUE)
data$ReportedEfficacyLin = scale(data$ReportedEfficacyLin, scale= FALSE, center = TRUE)
data$ReportedRewRateLin = scale(data$ReportedRewRateLin, scale= FALSE, center = TRUE)


# Z-score the efficacy and reward rate estimates within subject

# Previous trial
for (s in 1:length(unique(data$SubID))){
  data$mbased_efficacy_prev[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_efficacy_prev[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data$mbased_reward_prev[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_reward_prev[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

# Current trial
for (s in 1:length(unique(data$SubID))){
  data$mbased_efficacy[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_efficacy[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data$mbased_reward[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_reward[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data$mbased_efficacy_feedback_locked[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_efficacy_feedback_locked[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data$mbased_reward_feedback_locked[data$SubID == unique(data$SubID)[s]] = scale(data$mbased_reward_feedback_locked[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data$mbased_efficacy = data$mbased_efficacy/5.837704
data$mbased_efficacy_prev = data$mbased_efficacy_prev/5.837704
data$mbased_efficacy_feedback_locked = data$mbased_efficacy_feedback_locked/5.837704
data$mbased_reward = data$mbased_reward/7.99977
data$mbased_reward_prev = data$mbased_reward_prev/7.99977
data$mbased_reward_feedback_locked = data$mbased_reward_feedback_locked/7.99977


# Learning rates
data$LR_efficacy_raw = data$LR_efficacy
data$LR_reward_raw = data$LR_reward

data$LR_efficacy = scale(data$LR_efficacy, scale= FALSE, center = TRUE)
data$LR_reward = scale(data$LR_reward, scale= FALSE, center = TRUE)

data$signed_PE_efficacy = scale(data$signed_PE_efficacy, scale= FALSE, center = TRUE)
data$unsigned_PE_efficacy = scale(data$unsigned_PE_efficacy, scale= FALSE, center = TRUE)
data$signed_PE_reward = scale(data$signed_PE_reward, scale= FALSE, center = TRUE)
data$unsigned_PE_reward = scale(data$unsigned_PE_reward, scale= FALSE, center = TRUE)


# z-score the learning rates for efficacy
for (s in 1:length(unique(data$SubID))){
  data$LR_efficacy_raw_zscored[data$SubID == unique(data$SubID)[s]] = scale(data$LR_efficacy_raw[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = FALSE)
}

# z-score the learning rates for reward
for (s in 1:length(unique(data$SubID))){
  data$LR_reward_raw_zscored[data$SubID == unique(data$SubID)[s]] = scale(data$LR_reward_raw[data$SubID == unique(data$SubID)[s]], center = TRUE, scale = FALSE)
}


data$runAvgRewRate = scale(data$runAvgRewRate, scale= FALSE, center = TRUE)
data$runAvgEfficacy = scale(data$runAvgEfficacy, scale= FALSE, center = TRUE)

data$mbased_efficacy_update = scale(data$mbased_efficacy_update, scale= FALSE, center = TRUE)
data$mbased_efficacy_absolute_update = scale(data$mbased_efficacy_absolute_update, scale= FALSE, center = TRUE)
data$mbased_reward_update = scale(data$mbased_reward_prev, scale= FALSE, center = TRUE)
data$mbased_reward_absolute_update = scale(data$mbased_reward_absolute_update, scale= FALSE, center = TRUE)

data$Spearman_r_eff = scale(data$Spearman_r_eff, scale= FALSE, center = TRUE)
data$Spearman_r_rew_rate = scale(data$Spearman_r_rew_rate, scale= FALSE, center = TRUE)

# Contrast code the categorical predictors
data$feedbackOrder = as.factor(data$feedbackOrder)
contrasts(data$feedbackOrder) <- contr.sdif(2)
data$EffLvl = as.factor(data$EffLvl)
data$IsRewarded = as.factor(data$IsRewarded)
contrasts(data$EffLvl) <- contr.sdif(2)
colnames(attr(data$EffLvl, "contrasts")) =  c("NEff_min_Eff") # Change the name of the contrasts
contrasts(data$IsRewarded) <- contr.sdif(2)
colnames(attr(data$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data$Congruency = ordered(data$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data$Congruency) <- contr.sdif(3)
colnames(attr(data$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
data$Accuracy = data$Acc
data$Accuracy = ifelse(data$Accuracy==1,"Correct","Incorrect") 
data$Accuracy = as.factor(data$Accuracy)
contrasts(data$Accuracy) <- contr.sdif(2)
colnames(attr(data$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data$Rew_eff_fb = as.factor(data$Rew_eff_fb)
data$Rew_eff_fb = ordered(data$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data$Eff_rew_fb = as.factor(data$Eff_rew_fb)
data$Eff_rew_fb = ordered(data$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts

# contrasts(data$EffLvl) <- contr.sdif(2)
contrasts(data$Efficacy_T1) <- contr.sdif(2)
contrasts(data$Efficacy_T2) <- contr.sdif(2)
contrasts(data$Efficacy_T3) <- contr.sdif(2)
contrasts(data$Efficacy_T4) <- contr.sdif(2)
contrasts(data$Efficacy_T5) <- contr.sdif(2)
contrasts(data$Efficacy_T6) <- contr.sdif(2)
contrasts(data$Efficacy_T7) <- contr.sdif(2)
contrasts(data$Efficacy_T8) <- contr.sdif(2)
contrasts(data$Efficacy_T9) <- contr.sdif(2)
contrasts(data$Efficacy_T10) <- contr.sdif(2)

# contrasts(data$IsRewarded) <- contr.sdif(2)
contrasts(data$IsRewarded_T1) <- contr.sdif(2)
contrasts(data$IsRewarded_T2) <- contr.sdif(2)
contrasts(data$IsRewarded_T3) <- contr.sdif(2)
contrasts(data$IsRewarded_T4) <- contr.sdif(2)
contrasts(data$IsRewarded_T5) <- contr.sdif(2)
contrasts(data$IsRewarded_T6) <- contr.sdif(2)
contrasts(data$IsRewarded_T7) <- contr.sdif(2)
contrasts(data$IsRewarded_T8) <- contr.sdif(2)
contrasts(data$IsRewarded_T9) <- contr.sdif(2)
contrasts(data$IsRewarded_T10) <- contr.sdif(2)





#### LR & BIC correlation ####

library(Hmisc)

# Efficacy

# Positive learning rates
rcorr(learning_rates_efficacy$positive_learning_rate,learning_rates_efficacy$V6,type = "spearman")

PositiveLR<-ggplot(learning_rates_efficacy,aes(x= positive_learning_rate, y=V6)) +
  geom_point(size=3) +
  ggtitle('Efficacy learning rates and model fit')+
  ylab("BIC") +
  xlab("Positive learning rate") +
  # ylim(6.35,6.45)+ #edit this
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) 
PositiveLR

# Negative learning rates
rcorr(learning_rates_efficacy$negative_learning_rate,learning_rates_efficacy$V6,type = "spearman")

NegativeLR<-ggplot(learning_rates_efficacy,aes(x= negative_learning_rate, y=V6)) +
  geom_point(size=3) +
  ggtitle('Efficacy learning rates and model fit')+
  ylab("BIC") +
  xlab("Negative learning rate") +
  # ylim(6.35,6.45)+ #edit this
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) 
NegativeLR


# Reward

# Positive learning rates
rcorr(learning_rates_reward$positive_learning_rate,learning_rates_reward$V6,type = "spearman")

PositiveLR<-ggplot(learning_rates_reward,aes(x= positive_learning_rate, y=V6)) +
  geom_point(size=3) +
  ggtitle('Reward learning rates and model fit')+
  ylab("BIC") +
  xlab("Positive learning rate") +
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) 
PositiveLR

# Negative learning rates
rcorr(learning_rates_reward$negative_learning_rate,learning_rates_reward$V6,type = "spearman")

NegativeLR<-ggplot(learning_rates_reward,aes(x= negative_learning_rate, y=V6)) +
  geom_point(size=3) +
  ggtitle('Reward learning rates and model fit')+
  ylab("BIC") +
  xlab("Negative learning rate") +
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) 
NegativeLR













#### Efficacy function ####

  data.efficacy = data
  
  data.efficacy = data.efficacy %>%
    group_by(SubID) %>%
    mutate(quartile = ntile(mbased_efficacy_prev,4))
  
  data.efficacy$quartile = as.factor(data.efficacy$quartile) 
  
  df.quartiles = cbind(data.efficacy$SubID,data.efficacy$AccRT,data.efficacy$Acc,data.efficacy$quartile)
  colnames(df.quartiles) = c("SubID","AccRT","Acc","quartile")
  
  df.quartiles = as.data.frame(df.quartiles)
  df.quartiles = df.quartiles[complete.cases(df.quartiles[,4]), ]
  
  library(Rmisc)
  #accuracy and RT summary table
  summary <- Rmisc::summarySEwithin(df.quartiles, measurevar="AccRT",
                                    withinvars=c("quartile"),
                                    idvar="SubID", na.rm=TRUE, conf.interval=.95)
  
  summary = summary[-6,]
  
pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS8_A.pdf") 

  RT_Plot<-ggplot(summary,aes(x= quartile, y=AccRT)) +
    geom_line(size=1,position = position_jitter(width=0.1,seed = 123))+
    geom_point(size=3,position = position_jitter(width=0.1,seed = 123)) +
    geom_linerange(size=1, aes(ymin=AccRT-se,ymax=AccRT+se),position = position_jitter(width=0.1,seed = 123)) + 
    ggtitle('Reaction times')+
    ylab("Reaction times (ms)") +
    xlab("Model-based efficacy quartile") +
    ylim(640,660)+ #edit this
    theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          text = element_text(size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) +
    scale_color_manual(name="Penalty",values=c(rgb(0.1, 0.1, 0.1),rgb(0.5, 0.5, 0.5)))
  RT_Plot

  dev.off()
  RT_Plot
summary <- Rmisc::summarySEwithin(df.quartiles, measurevar="Acc",
                                  withinvars=c("quartile"),
                                  idvar="SubID", na.rm=TRUE, conf.interval=.95)

summary = summary[-6,]

pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS8_B.pdf") 


Acc_Plot<-ggplot(summary,aes(x= quartile, y=Acc)) +
  geom_line(size=1,position = position_jitter(width=0.1,seed = 123))+
  geom_point(size=3,position = position_jitter(width=0.1,seed = 123)) +
  geom_linerange(size=1, aes(ymin=Acc-se,ymax=Acc+se),position = position_jitter(width=0.1,seed = 123)) + 
  ggtitle('Accuracy')+
  ylab("Accuracy") +
  xlab("Model-based efficacy quartile") +
  scale_y_continuous(breaks=seq(0.75, 0.85, 0.05), limits=c(0.75, 0.85))+
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) +
  scale_color_manual(name="Penalty",values=c(rgb(0.1, 0.1, 0.1),rgb(0.5, 0.5, 0.5)))
Acc_Plot
dev.off()
Acc_Plot















#### CNV raw data ####

data.efficacy = data=data[!is.na(data$Resp) & data$RT>200 & !is.na(data$mbased_reward_prev),]

data.efficacy = data.efficacy %>%
  group_by(SubID) %>%
  mutate(quartile = ntile(mbased_efficacy_prev,4))

data.efficacy$quartile = as.factor(data.efficacy$quartile) 

df.quartiles = cbind(data.efficacy$SubID,data.efficacy$CNV10001500,data.efficacy$Acc,data.efficacy$quartile)
colnames(df.quartiles) = c("SubID","CNV","Acc","quartile")

df.quartiles = as.data.frame(df.quartiles)
df.quartiles = df.quartiles[complete.cases(df.quartiles[,4]), ]

library(Rmisc)
#accuracy and RT summary table
summary <- Rmisc::summarySEwithin(df.quartiles, measurevar="CNV",
                                  withinvars=c("quartile"),
                                  idvar="SubID", na.rm=TRUE, conf.interval=.95)

summary = summary[-6,]

pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS8_E.pdf") 

CNV_Plot<-ggplot(summary,aes(x= quartile, y=CNV)) +
  geom_line(size=1,position = position_jitter(width=0.1,seed = 123))+
  geom_point(size=3,position = position_jitter(width=0.1,seed = 123)) +
  geom_linerange(size=1, aes(ymin=CNV-se,ymax=CNV+se),position = position_jitter(width=0.1,seed = 123)) + 
  ggtitle('CNV')+
  ylab("CNV") +
  xlab("Model-based efficacy quartile") +
  # ylim(6.35,6.45)+ #edit this
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) +
  scale_color_manual(name="Penalty",values=c(rgb(0.1, 0.1, 0.1),rgb(0.5, 0.5, 0.5)))
CNV_Plot

dev.off()
CNV_Plot



#### P3b raw data ####

data.efficacy = data[data$Acc == 1 & data$RT>200,]

data.efficacy = data.efficacy %>%
  group_by(SubID) %>%
  mutate(quartile = ntile(mbased_efficacy_prev,4))

data.efficacy$quartile = as.factor(data.efficacy$quartile) 

df.quartiles = cbind(data.efficacy$SubID,data.efficacy$P3bE350500,data.efficacy$Acc,data.efficacy$quartile,data.efficacy$EffLvl)
colnames(df.quartiles) = c("SubID","P3b","Acc","quartile","EfficacyFeedback")

df.quartiles = as.data.frame(df.quartiles)
df.quartiles = df.quartiles[complete.cases(df.quartiles[,4]), ]

library(Rmisc)
#accuracy and RT summary table
summary <- Rmisc::summarySEwithin(df.quartiles, measurevar="P3b",
                                  withinvars=c("quartile","EfficacyFeedback"),
                                  idvar="SubID", na.rm=TRUE, conf.interval=.95)

# summary = summary[-6,]

pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS8_F.pdf") 

P3b_Plot<-ggplot(summary,aes(x= quartile, y=P3b,color=EfficacyFeedback)) +
  geom_line(size=1,position = position_jitter(width=0.1,seed = 123))+
  geom_point(size=3,position = position_jitter(width=0.1,seed = 123)) +
  geom_linerange(size=1, aes(ymin=P3b-se,ymax=P3b+se),position = position_jitter(width=0.1,seed = 123)) + 
  ggtitle('P3b')+
  ylab("P3b") +
  xlab("Model-based efficacy quartile") +
  # ylim(6.35,6.45)+ #edit this
  theme(axis.line = element_line(colour = "black"),plot.margin=grid::unit(c(5,15,5,5), "mm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 20, face = "bold")) +
  scale_color_manual(name="Efficacy feedback",values=c(rgb(0.1, 0.1, 0.1),rgb(0.5, 0.5, 0.5)))
P3b_Plot

dev.off()
P3b_Plot



#### brms settings ####


#help stan run faster
# rstan_options(auto_write = TRUE)
cores=options(mc.cores = parallel::detectCores())
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')


#### Reward probability yoking ####
# data$Acc = as.numeric(data$Acc)
# # Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Model efficacy binary
model.reward_probability_EffLvl = brm(IsRewarded ~ EffLvl + (EffLvl|SubID),
                                         data=data,
                                         family="bernoulli",
                                         iter  = 10000,
                                         warmup = 8000,
                                         save_pars=save_pars(all=TRUE),
                                         control = list(adapt_delta = 0.99),
                                         sample_prior = TRUE)
saveRDS(model.reward_probability_EffLvl,file="model.model.reward_probability_EffLvl1.rds")

# Model efficacy model-based
model.reward_probability_Effmbased = brm(IsRewarded ~ mbased_efficacy_feedback_locked + (mbased_efficacy_feedback_locked|SubID),
                                      data=data,
                                      family="bernoulli",
                                      iter  = 10000,
                                      warmup = 8000,
                                      save_pars=save_pars(all=TRUE),
                                      control = list(adapt_delta = 0.99),
                                      sample_prior = TRUE)
saveRDS(model.reward_probability_Effmbased,file="model.model.reward_probability_Effmbased1.rds")

# Model efficacy model-based predicts reward model-based
model.Rewmbased_Effmbased = brm(mbased_reward ~ mbased_efficacy + (mbased_efficacy|SubID),
                                         data=data,
                                         family=gaussian(),
                                         iter  = 10000,
                                         warmup = 8000,
                                         save_pars=save_pars(all=TRUE),
                                         control = list(adapt_delta = 0.99),
                                         sample_prior = TRUE)
saveRDS(model.Rewmbased_Effmbased,file="model.model.Rewmbased_Effmbased1.rds")


#### Subjective reward rate ####

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors
prior = c(
  prior(normal(0.5,0.2), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 0.2), class = b)) # a wide prior

# Predicting reward probe response by reward and efficacy feedback
reward = brm(RewRateProbeResp ~  mbased_reward + mbased_efficacy + 
                 (mbased_reward + mbased_efficacy|SubID),
               data=data[!is.na(data$RewRateProbeResp),],
               family=gaussian(),
               prior = prior,
               iter  = 20000,
               warmup = 19000,
               save_pars=save_pars(all=TRUE),
               control = list(adapt_delta = 0.99),
               sample_prior = TRUE)
saveRDS(reward,file="model.reward.subest_by_modelbased.rds")

# # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(0.5,0.2), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 0.2), class = b)) # a wide prior
# 
# # Predicting reward probe response by reward and efficacy feedback
# reward = brm(RewRateProbeResp ~  1 + IsRewarded +
#                IsRewarded_T1 +
#                IsRewarded_T2 +
#                IsRewarded_T3 +
#                IsRewarded_T4 +
#                EffLvl +
#                Efficacy_T1 +
#                Efficacy_T2 +
#                Efficacy_T3 +
#                Efficacy_T4 +
#                  (1|SubID),
#                data=data[!is.na(data$RewRateProbeResp),],
#                family=gaussian(),
#                prior = prior,
#                iter  = 20000,
#                warmup = 19000,
#                save_pars=save_pars(all=TRUE),
#                control = list(adapt_delta = 0.99),
#                sample_prior = TRUE)
# saveRDS(reward,file="model.reward.previous.efficacy_5.rds")


#### Subjective efficacy ####

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors
prior = c(
  prior(normal(0.5,0.2), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 0.2), class = b)) # a wide prior

# Predicting efficacy probe response by reward and efficacy feedback
reward = brm(EfficacyProbeResp ~  mbased_reward + mbased_efficacy + 
               (mbased_reward + mbased_efficacy|SubID),
             data=data[!is.na(data$EfficacyProbeResp),],
             family=gaussian(),
             prior = prior,
             iter  = 20000,
             warmup = 19000,
             save_pars=save_pars(all=TRUE),
             control = list(adapt_delta = 0.99),
             sample_prior = TRUE)
saveRDS(reward,file="model.efficacy.subest_by_modelbased.rds")

# # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(0.5,0.2), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 0.2), class = b)) # a wide prior
# 
# # Predicting efficacy probe response by efficacy and reward feedback
# efficacy = brm(EfficacyProbeResp ~  1 + 
#                  EffLvl +
#                  Efficacy_T1 +
#                  Efficacy_T2 +
#                  Efficacy_T3 +
#                  Efficacy_T4 +
#                  IsRewarded +
#                  IsRewarded_T1 +
#                  IsRewarded_T2 +
#                  IsRewarded_T3 +
#                  IsRewarded_T4 +
#                  (1|SubID),
#                data=data[!is.na(data$EfficacyProbeResp),],
#                family=gaussian(),
#                prior = prior,
#                iter  = 20000,
#                warmup = 19000,
#                save_pars=save_pars(all=TRUE),
#                control = list(adapt_delta = 0.99),
#                sample_prior = TRUE)
# saveRDS(efficacy,file="model.efficacy.by.previous.reward_5.rds")



#### t-test for LRs ####

# ## Efficacy
# 
# # wide to long format
# lr_ttest = subset(learning_rates_efficacy, select=c("positive_learning_rate","negative_learning_rate","SubID"))
# lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)
# 
# # set the variable types
# lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
# contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)
# 
# # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(0.5,0.5), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 0.5), class = b)) # a wide prior
# 
# # All subjects
# m = brm(
#   Value ~ Learning_rate_type, sigma ~ Learning_rate_type,
#   prior=prior,
#   family=student,
#   data = lr_ttest,
#   iter  = 20000,
#   warmup = 19000,
#   save_pars=save_pars(all=TRUE),
#   control = list(adapt_delta = 0.99),
#   sample_prior = TRUE)
# saveRDS(m,file="model.ttest.learning_rates_efficacy.rds")
# 
# 
# # Even subjects
# 
# # wide to long format
# lr_ttest = subset(learning_rates_efficacy, select=c("positive_learning_rate","negative_learning_rate","SubID"))
# # take only even subjects
# lr_ttest = subset(learning_rates_efficacy, SubID%% 2 ==0)
# 
# lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)
# 
# # set the variable types
# lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
# contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)
# 
# # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(0.5,0.5), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 0.5), class = b)) # a wide prior
# 
# m = brm(
#   Value ~ Learning_rate_type, sigma ~ Learning_rate_type,
#   prior=prior,
#   family=student,
#   data = lr_ttest,
#   iter  = 20000,
#   warmup = 19000,
#   save_pars=save_pars(all=TRUE),
#   control = list(adapt_delta = 0.99),
#   #cores = 12,
#   sample_prior = TRUE)
# saveRDS(m,file="model.ttest.learning_rates_efficacy_even_subs.rds")


# Odd subjects

# wide to long format
lr_ttest = subset(learning_rates_efficacy, select=c("positive_learning_rate","negative_learning_rate","SubID"))
# take only odd subjects
lr_ttest = subset(learning_rates_efficacy, SubID%% 2 !=0)

lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)

# set the variable types
lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

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
  save_pars=save_pars(all=TRUE),
  sample_prior = TRUE)
saveRDS(m,file="model.ttest.learning_rates_efficacy_odd_subs.rds")


## Reward

# wide to long format
lr_ttest = subset(learning_rates_reward, select=c("positive_learning_rate","negative_learning_rate","SubID"))
lr_ttest = lr_ttest %>% gather(Learning_rate_type,Value,positive_learning_rate,negative_learning_rate)

# set the variable types
lr_ttest$Learning_rate_type = as.factor(lr_ttest$Learning_rate_type)
contrasts(lr_ttest$Learning_rate_type) <- contr.sdif(2)

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

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
  save_pars=save_pars(all=TRUE),
  sample_prior = TRUE)
saveRDS(m,file="model.ttest.learning_rates_reward.rds")




#### Reaction times ####

# # Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors (Based on Fromero et al., 2020)
prior = c(
  prior(normal(624, 58.69), class = Intercept),
  prior(normal(15.54, 21.93), class = b, coef=CongruencyFacilitation),
  prior(normal(61.12, 37.49), class = b, coef=CongruencyInterference),
  prior(normal(-10.26, 19.51), class = b, coef=mbased_efficacy_prev),
  prior(normal(0, 19.51), class = b, coef=mbased_reward_prev))

# Congruency + Efficacy model
model.congruency_plus_efficacy.RT = brm(AccRT ~ Congruency + mbased_efficacy_prev + mbased_reward_prev + (Congruency + mbased_efficacy_prev + mbased_reward_prev|SubID),
                          data=data[!is.na(data$mbased_reward_prev) & !is.na(data$AccRT),],
                          family=exgaussian(),
                          prior = prior,
                          iter  = 20000,
                          warmup = 18000,
                          save_pars=save_pars(all=TRUE),
                          control = list(adapt_delta = 0.99),
                          sample_prior = TRUE)
saveRDS(model.congruency_plus_efficacy.RT,file="model.congruency_plus_efficacy.RT.rds")


#### Accuracy ####
data$Acc = as.numeric(data$Acc)
# # Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors (Based on Fromero et al., 2020)
prior = c(
  prior(normal(2.11, 0.81), class = Intercept),
  prior(normal(-0.45, 0.64), class = b, coef=CongruencyFacilitation),
  prior(normal(-0.53, 0.81), class = b, coef=CongruencyInterference),
  prior(normal(0.09, 0.32), class = b, coef=mbased_efficacy_prev),
  prior(normal(0, 0.32), class = b, coef=mbased_reward_prev))

# Congruency + Efficacy model
model.congruency_plus_efficacy.Acc = brm(Acc ~ Congruency + mbased_efficacy_prev + mbased_reward_prev + (Congruency + mbased_efficacy_prev + mbased_reward_prev|SubID),
                                        data=data[!is.na(data$Resp) & !is.na(data$mbased_reward_prev),],
                                        family="bernoulli",
                                        prior = prior,
                                        iter  = 20000,
                                        warmup = 18000,
                                        save_pars=save_pars(all=TRUE),
                                        control = list(adapt_delta = 0.99),
                                        sample_prior = TRUE)
saveRDS(model.congruency_plus_efficacy.Acc,file="model.congruency_plus_efficacy.Acc.rds")

#### Late CNV ####

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors (Based on Fromero et al., 2020 )
prior = c(
  prior(normal(-0.16, 1.59), class = Intercept),
  prior(normal(-0.30, 0.73), class = b, coef=mbased_efficacy_prev),
  prior(normal(0, 0.73), class = b, coef=mbased_reward_prev))

# Efficacy estimate model
model.efficacy.lateCNV = brm(CNV10001500 ~ mbased_efficacy_prev + mbased_reward_prev + (mbased_efficacy_prev + mbased_reward_prev|SubID),
                             data=data[!is.na(data$Resp) & data$RT>200 & !is.na(data$mbased_reward_prev),],
                             family=gaussian(),
                             prior = prior,
                             iter  = 20000,
                             warmup = 19000,
                             save_pars=save_pars(all=TRUE),
                             control = list(adapt_delta = 0.99),
                             sample_prior = TRUE)
saveRDS(model.efficacy.lateCNV,file="model.efficacy.lateCNV.rds")



#### Reaction times predicted by the late CNV ####

# # # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(650, 200), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 50), class = b)) # a wide prior
# 
# data.CNV = data
# 
# # z-scoreand center the CNV
# for (s in 1:length(unique(data$SubID))){
#   data.CNV$CNV10001500[data$SubID == unique(data$SubID)[s]] = scale(data.CNV$CNV10001500[data.CNV$SubID == unique(data.CNV$SubID)[s]], center = TRUE, scale = TRUE)
# }
# 
# # The model
# model.RT.lateCNV = brm(AccRT ~ CNV10001500 + (CNV10001500|SubID),
#                        data=data.CNV[!is.na(data.CNV$AccRT),],
#                        family=exgaussian(),
#                        prior = prior,
#                        iter  = 20000,
#                        warmup = 18000,
#                        save_pars=save_pars(all=TRUE),
#                        control = list(adapt_delta = 0.99),
#                        sample_prior = TRUE)
# saveRDS(model.RT.lateCNV,file="model.RT.lateCNV.rds")
# 
# 
# model.RT.lateCNV = readRDS("model.RT.lateCNV.rds")


#### Accuracy predicted by the late CNV####

# # # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# # Set the priors
# prior = c(
#   prior(normal(0.7, 0.2), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 0.2), class = b)) # a wide prior
# 
# # The model
# model.Acc.lateCNV = brm(Acc ~ CNV10001500 + (CNV10001500|SubID),
#                                           data=data.CNV[!is.na(data.CNV$Resp),],
#                                           family="bernoulli",
#                                           prior = prior,
#                                           iter  = 20000,
#                                           warmup = 18000,
#                                           save_pars=save_pars(all=TRUE),
#                                           control = list(adapt_delta = 0.99),
#                                           sample_prior = TRUE)
# saveRDS(model.Acc.lateCNV,file="model.Acc.lateCNV.rds")


#### P3b - Reward ####

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors
prior = c(
  prior(normal(5, 5), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 5), class = b)) # a wide prior


# Estimate times feedback model
model.reward.P3b.reward_times_feedback = brm(P3bR350500 ~ mbased_reward_prev * IsRewarded  + mbased_efficacy_feedback_locked + (mbased_reward_prev * IsRewarded  + mbased_efficacy_feedback_locked|SubID),
                                             data=data[data$Acc == 1 & data$RT>200,],
                                             family=gaussian(),
                                             prior = prior,
                                             iter  = 20000,
                                             warmup = 19000,
                                             save_pars=save_pars(all=TRUE),
                                             control = list(adapt_delta = 0.99),
                                             sample_prior = TRUE)
saveRDS(model.reward.P3b.reward_times_feedback,file="model.reward.P3b.reward_times_feedback.rds")


# Unsigned PE times reward feedback model
model.reward.P3b.PE_times_LR = brm(P3bR350500 ~ unsigned_PE_reward * LR_reward_raw_zscored  + mbased_efficacy_feedback_locked + (unsigned_PE_reward * LR_reward_raw_zscored  + mbased_efficacy_feedback_locked|SubID),
                                   data=data[data$Acc == 1 & data$RT>200,],
                                   family=gaussian(),
                                   prior = prior,
                                   iter  = 20000,
                                   warmup = 19000,
                                   save_pars=save_pars(all=TRUE),
                                   control = list(adapt_delta = 0.99),
                                   sample_prior = TRUE)
saveRDS(model.reward.P3b.PE_times_LR,file="model.reward.P3b.PE_times_LR.rds")


#### P3b - Efficacy ####

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))

# Set the priors
prior = c(
  prior(normal(5, 3), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 3), class = b)) # a wide prior

# Unsigned PE times efficacy feedback model
model.efficacy.P3b.PE_times_LR = brm(P3bE350500 ~ unsigned_PE_efficacy * LR_efficacy_raw_zscored  + mbased_reward_feedback_locked + (unsigned_PE_efficacy * LR_efficacy_raw_zscored  + mbased_reward_feedback_locked|SubID),
                                     data=data[data$Acc == 1 & data$RT>200,],
                                     family=gaussian(),
                                     prior = prior,
                                     iter  = 20000,
                                     warmup = 19000,
                                     save_pars=save_pars(all=TRUE),
                                     control = list(adapt_delta = 0.99),
                                     sample_prior = TRUE)
saveRDS(model.efficacy.P3b.PE_times_LR,file="model.efficacy.P3b.PE_times_LR.rds")

# Estimate times feedback model
model.efficacy.P3b.efficacy_times_feedback = brm(P3bE350500 ~ mbased_efficacy_prev * EffLvl  + mbased_reward_feedback_locked + (mbased_efficacy_prev * EffLvl  + mbased_reward_feedback_locked|SubID),
                                                 data=data[data$Acc == 1 & data$RT>200,],
                                                 family=gaussian(),
                                                 prior = prior,
                                                 iter  = 20000,
                                                 warmup = 19000,
                                                 save_pars=save_pars(all=TRUE),
                                                 control = list(adapt_delta = 0.99),
                                                 sample_prior = TRUE)
saveRDS(model.efficacy.P3b.efficacy_times_feedback,file="model.efficacy.P3b.efficacy_times_feedback.rds")

#### Topographies ####


#### P3b efficacy  - estimate * feedback ####

# Topography
data_topo  =  cbind(data.raw, P3bE350500) # the size needs to be exactly the same as the raw eeg data

# delete RTs below 200ms
data_topo$RT = ifelse(data_topo$RT<200,NA,data_topo$RT)

# Rescale the data_topo for the regressions
# Previous trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Current trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data_topo$mbased_efficacy = data_topo$mbased_efficacy/5.837704
data_topo$mbased_efficacy_prev = data_topo$mbased_efficacy_prev/5.837704
data_topo$mbased_efficacy_feedback_locked = data_topo$mbased_efficacy_feedback_locked/5.837704
data_topo$mbased_reward = data_topo$mbased_reward/7.99977
data_topo$mbased_reward_prev = data_topo$mbased_reward_prev/7.99977
data_topo$mbased_reward_feedback_locked = data_topo$mbased_reward_feedback_locked/7.99977


data_topo$feedbackOrder = ifelse(data_topo$feedbackOrder==1,"Efficacy_First","Efficacy_Second") # 1 is first efficacy
data_topo$feedbackOrder = as.factor(data_topo$feedbackOrder)
contrasts(data_topo$feedbackOrder) <- contr.sdif(2)
data_topo$EffLvl = as.factor(data_topo$EffLvl)
data_topo$IsRewarded = as.factor(data_topo$IsRewarded)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
colnames(attr(data_topo$EffLvl, "contrasts")) =  c("Eff_min_NEff") # Change the name of the contrasts
contrasts(data_topo$IsRewarded) <- contr.sdif(2)
colnames(attr(data_topo$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data_topo$Congruency = ordered(data_topo$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data_topo$Congruency) <- contr.sdif(3)
colnames(attr(data_topo$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
# data_topo$Accuracy = data_topo$Acc
# data_topo$Accuracy = ifelse(data_topo$Accuracy==1,"Correct","Incorrect") 
# data_topo$Accuracy = as.factor(data_topo$Accuracy)
# contrasts(data_topo$Accuracy) <- contr.sdif(2)
colnames(attr(data_topo$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data_topo$Rew_eff_fb = as.factor(data_topo$Rew_eff_fb)
data_topo$Rew_eff_fb = ordered(data_topo$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data_topo$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data_topo$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data_topo$Eff_rew_fb = as.factor(data_topo$Eff_rew_fb)
data_topo$Eff_rew_fb = ordered(data_topo$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data_topo$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data_topo$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts
data_topo$LR = scale(data_topo$LR, scale= FALSE, center = TRUE)

# estimates
plotEmat = as.data.frame(matrix(nrow = 65, ncol=5))  #as many rows as the number of the electrodes, ncol is the number of fixed effects in the model

# Set the priors 
prior = c(
  prior(normal(5, 3), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 3), class = b)) # a wide prior

# # Take only the correct and not too fast responses
# data_topo = data_topo[data_topo$Acc == 1 & data_topo$RT>200,]
# 
# # Create a list of data
# data_topo_list = list()
# 
# # Set the working directory where to save the models
# setwd(here("Analyses/Stats/brms/brms_models"))
# 
# j=1
# for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependent variable
#   data_topo$DV  <- data_topo[,i]
#   data_topo_list[[j]] = data_topo
#   j=j+1
# }
# 
# # Fit a model to each of the datasets
# Plotmod <- 
#   brm_multiple(DV ~ EffLvl * mbased_efficacy_prev + mbased_reward_feedback_locked + (1|SubID),
#       data=data_topo_list,
#       family=gaussian(),
#       prior = prior,
#       iter  = 30,
#       warmup=20,
#       combine = FALSE,
#       )
# 
# saveRDS(Plotmod,file="Test_Topo.rds")



# fit a model to each electrode and extract the estimates 

# Fit the first model
data_topo$DV  <- data_topo[,ncol(data.raw)+1]
Plotmod <- 
  brm(DV ~ EffLvl * mbased_efficacy_prev + mbased_reward_feedback_locked + (1|SubID),
      data=data_topo[data_topo$Acc == 1 & data_topo$RT>200,],
      family=gaussian(),
      prior = prior,
      iter  = 3000,
      warmup=2000)

k <- summary(Plotmod)
f <- dimnames(k$fixed)[[1]]
colnames(plotEmat) <-  f
plotEmat[1,] <- k$fixed[,1] #estimate
# Fit the rest of the models
j=2
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependent variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod1 <- 
    update(Plotmod,newdata=data_topo[data_topo$Acc == 1 & data_topo$RT>200,])
  
  k <- summary(Plotmod1)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
  # to prevent crashes
  Sys.sleep(30)
  rm(Plotmod1)
  # gc()
}

# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_reward_feedback_locked", "Rew_estimate")
colnames(plotEmat) <- sedit(colnames(plotEmat), "EffLvlEff_min_NEffbymbased_efficacy_prev", "Eff_min_NEffbyEff_estimate")

# Save as .mat files
writeMat('C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Data/EEG/Export/plotesP3beff.mat',EmatP3beff=plotEmat)


#### P3b efficacy  - PE * LR ####

# Topography
data_topo  =  cbind(data.raw, P3bE350500) # the size needs to be exactly the same as the raw eeg data

# delete RTs below 200ms
data_topo$RT = ifelse(data_topo$RT<200,NA,data_topo$RT)

data_topo$mbased_reward_feedback_locked = scale(data_topo$mbased_reward_feedback_locked, scale= FALSE, center = TRUE)

# Rescale the data_topo for the regressions
# Previous trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Current trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data_topo$mbased_efficacy = data_topo$mbased_efficacy/5.837704
data_topo$mbased_efficacy_prev = data_topo$mbased_efficacy_prev/5.837704
data_topo$mbased_efficacy_feedback_locked = data_topo$mbased_efficacy_feedback_locked/5.837704
data_topo$mbased_reward = data_topo$mbased_reward/7.99977
data_topo$mbased_reward_prev = data_topo$mbased_reward_prev/7.99977
data_topo$mbased_reward_feedback_locked = data_topo$mbased_reward_feedback_locked/7.99977

data_topo$feedbackOrder = ifelse(data_topo$feedbackOrder==1,"Efficacy_First","Efficacy_Second") # 1 is first efficacy
data_topo$feedbackOrder = as.factor(data_topo$feedbackOrder)
contrasts(data_topo$feedbackOrder) <- contr.sdif(2)
data_topo$EffLvl = as.factor(data_topo$EffLvl)
data_topo$IsRewarded = as.factor(data_topo$IsRewarded)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
colnames(attr(data_topo$EffLvl, "contrasts")) =  c("Eff_min_NEff") # Change the name of the contrasts
contrasts(data_topo$IsRewarded) <- contr.sdif(2)
colnames(attr(data_topo$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data_topo$Congruency = ordered(data_topo$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data_topo$Congruency) <- contr.sdif(3)
colnames(attr(data_topo$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
data_topo$Accuracy = data_topo$Acc
data_topo$Accuracy = ifelse(data_topo$Accuracy==1,"Correct","Incorrect") 
data_topo$Accuracy = as.factor(data_topo$Accuracy)
contrasts(data_topo$Accuracy) <- contr.sdif(2)
colnames(attr(data_topo$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data_topo$Rew_eff_fb = as.factor(data_topo$Rew_eff_fb)
data_topo$Rew_eff_fb = ordered(data_topo$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data_topo$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data_topo$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data_topo$Eff_rew_fb = as.factor(data_topo$Eff_rew_fb)
data_topo$Eff_rew_fb = ordered(data_topo$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data_topo$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data_topo$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts
data_topo$LR = scale(data_topo$LR, scale= FALSE, center = TRUE)
data_topo$signed_PE_efficacy = scale(data_topo$signed_PE_efficacy, scale= FALSE, center = TRUE)
data_topo$unsigned_PE_efficacy = scale(data_topo$unsigned_PE_efficacy, scale= FALSE, center = TRUE)
data_topo$signed_PE_reward = scale(data_topo$signed_PE_reward, scale= FALSE, center = TRUE)
data_topo$unsigned_PE_reward = scale(data_topo$unsigned_PE_reward, scale= FALSE, center = TRUE)

data_topo$LR_efficacy_raw = data_topo$LR_efficacy

# z-score the learning rates for efficacy
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$LR_efficacy_raw_zscored[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$LR_efficacy_raw[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = FALSE)
}

# estimates
plotEmat = as.data.frame(matrix(nrow = 65, ncol=5))  #as many rows as the number of the electrodes, ncol is the number of fixed effects in the model

# Set the priors 
prior = c(
  prior(normal(5, 3), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 3), class = b)) # a wide prior

# fit a model to each electrode and extract the estimates 
j=1
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependent variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod <- 
    brm(DV ~ unsigned_PE_efficacy * LR_efficacy_raw_zscored  + mbased_reward_feedback_locked + (1|SubID),
        data=data_topo[data_topo$Acc == 1 & data_topo$RT>200,],
        family=gaussian(),
        prior = prior,
        iter  = 3000,
        warmup=2000)
  
  k <- summary(Plotmod)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
}

# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_reward_feedback_locked", "Rew_estimate")
colnames(plotEmat) <- sedit(colnames(plotEmat), "unsigned_PE_efficacy", "PE")
colnames(plotEmat) <- sedit(colnames(plotEmat), "LR_efficacy_raw_zscored", "LR")


# Save as .mat files
writeMat('C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Data/EEG/Export/plotesP3beff_PE_LR.mat',EmatP3beff=plotEmat)



#### P3b reward  - estimate * feedback ####

# Topography
data_topo  =  cbind(data.raw, P3bR350500) # the size needs to be exactly the same as the raw eeg data

# delete RTs below 200ms
data_topo$RT = ifelse(data_topo$RT<200,NA,data_topo$RT)

data_topo$mbased_efficacy_feedback_locked = scale(data_topo$mbased_efficacy_feedback_locked, scale= FALSE, center = TRUE)


# Rescale the data_topo for the regressions
# Previous trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Current trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data_topo$mbased_efficacy = data_topo$mbased_efficacy/5.837704
data_topo$mbased_efficacy_prev = data_topo$mbased_efficacy_prev/5.837704
data_topo$mbased_efficacy_feedback_locked = data_topo$mbased_efficacy_feedback_locked/5.837704
data_topo$mbased_reward = data_topo$mbased_reward/7.99977
data_topo$mbased_reward_prev = data_topo$mbased_reward_prev/7.99977
data_topo$mbased_reward_feedback_locked = data_topo$mbased_reward_feedback_locked/7.99977


data_topo$feedbackOrder = ifelse(data_topo$feedbackOrder==1,"Efficacy_First","Efficacy_Second") # 1 is first efficacy
data_topo$feedbackOrder = as.factor(data_topo$feedbackOrder)
contrasts(data_topo$feedbackOrder) <- contr.sdif(2)
data_topo$EffLvl = as.factor(data_topo$EffLvl)
data_topo$IsRewarded = as.factor(data_topo$IsRewarded)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
colnames(attr(data_topo$EffLvl, "contrasts")) =  c("Eff_min_NEff") # Change the name of the contrasts
contrasts(data_topo$IsRewarded) <- contr.sdif(2)
colnames(attr(data_topo$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data_topo$Congruency = ordered(data_topo$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data_topo$Congruency) <- contr.sdif(3)
colnames(attr(data_topo$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
data_topo$Accuracy = data_topo$Acc
data_topo$Accuracy = ifelse(data_topo$Accuracy==1,"Correct","Incorrect") 
data_topo$Accuracy = as.factor(data_topo$Accuracy)
contrasts(data_topo$Accuracy) <- contr.sdif(2)
colnames(attr(data_topo$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data_topo$Rew_eff_fb = as.factor(data_topo$Rew_eff_fb)
data_topo$Rew_eff_fb = ordered(data_topo$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data_topo$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data_topo$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data_topo$Eff_rew_fb = as.factor(data_topo$Eff_rew_fb)
data_topo$Eff_rew_fb = ordered(data_topo$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data_topo$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data_topo$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts
data_topo$LR = scale(data_topo$LR, scale= FALSE, center = TRUE)

# estimates
plotEmat = as.data.frame(matrix(nrow = 65, ncol=5))  #as many rows as the number of the electrodes, ncol is the number of fixed effects in the model

# Set the priors 
prior = c(
  prior(normal(5, 3), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 3), class = b)) # a wide prior

# fit a model to each electrode and extract the estimates 
j=1
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependet variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod <- 
    brm(DV ~ mbased_reward_prev * IsRewarded  + mbased_efficacy_feedback_locked + (1|SubID),
        data=data_topo[data_topo$Acc == 1 & data_topo$RT>200,],
        family=gaussian(),
        prior = prior,
        iter  = 3000,
        warmup=2000)
  
  k <- summary(Plotmod)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
}

# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_efficacy_feedback_locked", "Eff_estimate")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_reward_prevbyIsRewardedRew_min_NRew", "Rew_min_NRewbyRew_estimate")

# Save as .mat files
writeMat('C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Data/EEG/Export/plotesP3brew.mat',EmatP3brew=plotEmat)


#### P3b reward  - PE * LR ####

# Topography
data_topo  =  cbind(data.raw, P3bR350500) # the size needs to be exactly the same as the raw eeg data

# delete RTs below 200ms
data_topo$RT = ifelse(data_topo$RT<200,NA,data_topo$RT)

data_topo$mbased_efficacy_feedback_locked = scale(data_topo$mbased_efficacy_feedback_locked, scale= FALSE, center = TRUE)

# Rescale the data_topo for the regressions
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Current trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data_topo$mbased_efficacy = data_topo$mbased_efficacy/5.837704
data_topo$mbased_efficacy_prev = data_topo$mbased_efficacy_prev/5.837704
data_topo$mbased_efficacy_feedback_locked = data_topo$mbased_efficacy_feedback_locked/5.837704
data_topo$mbased_reward = data_topo$mbased_reward/7.99977
data_topo$mbased_reward_prev = data_topo$mbased_reward_prev/7.99977
data_topo$mbased_reward_feedback_locked = data_topo$mbased_reward_feedback_locked/7.99977

data_topo$feedbackOrder = ifelse(data_topo$feedbackOrder==1,"Efficacy_First","Efficacy_Second") # 1 is first efficacy
data_topo$feedbackOrder = as.factor(data_topo$feedbackOrder)
contrasts(data_topo$feedbackOrder) <- contr.sdif(2)
data_topo$EffLvl = as.factor(data_topo$EffLvl)
data_topo$IsRewarded = as.factor(data_topo$IsRewarded)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
colnames(attr(data_topo$EffLvl, "contrasts")) =  c("Eff_min_NEff") # Change the name of the contrasts
contrasts(data_topo$IsRewarded) <- contr.sdif(2)
colnames(attr(data_topo$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data_topo$Congruency = ordered(data_topo$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data_topo$Congruency) <- contr.sdif(3)
colnames(attr(data_topo$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
data_topo$Accuracy = data_topo$Acc
data_topo$Accuracy = ifelse(data_topo$Accuracy==1,"Correct","Incorrect") 
data_topo$Accuracy = as.factor(data_topo$Accuracy)
contrasts(data_topo$Accuracy) <- contr.sdif(2)
colnames(attr(data_topo$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data_topo$Rew_eff_fb = as.factor(data_topo$Rew_eff_fb)
data_topo$Rew_eff_fb = ordered(data_topo$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data_topo$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data_topo$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data_topo$Eff_rew_fb = as.factor(data_topo$Eff_rew_fb)
data_topo$Eff_rew_fb = ordered(data_topo$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data_topo$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data_topo$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts
data_topo$LR = scale(data_topo$LR, scale= FALSE, center = TRUE)

data_topo$signed_PE_efficacy = scale(data_topo$signed_PE_efficacy, scale= FALSE, center = TRUE)
data_topo$unsigned_PE_efficacy = scale(data_topo$unsigned_PE_efficacy, scale= FALSE, center = TRUE)
data_topo$signed_PE_reward = scale(data_topo$signed_PE_reward, scale= FALSE, center = TRUE)
data_topo$unsigned_PE_reward = scale(data_topo$unsigned_PE_reward, scale= FALSE, center = TRUE)

data_topo$LR_reward_raw = data_topo$LR_reward

# z-score the learning rates for efficacy
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$LR_reward_raw_zscored[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$LR_reward_raw[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = FALSE)
}

# estimates
plotEmat = as.data.frame(matrix(nrow = 65, ncol=5))  #as many rows as the number of the electrodes, ncol is the number of fixed effects in the model

# Set the priors 
prior = c(
  prior(normal(5, 3), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 3), class = b)) # a wide prior

# fit a model to each electrode and extract the estimates 
j=1
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependent variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod <- 
    brm(DV ~ unsigned_PE_reward * LR_reward_raw_zscored  + mbased_efficacy_feedback_locked + (1|SubID),
        data=data_topo[data_topo$Acc == 1 & data_topo$RT>200,],
        family=gaussian(),
        prior = prior,
        iter  = 3000,
        warmup=2000)
  
  k <- summary(Plotmod)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
}

# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_efficacy_feedback_locked", "Eff_estimate")
colnames(plotEmat) <- sedit(colnames(plotEmat), "unsigned_PE_reward", "PE")
colnames(plotEmat) <- sedit(colnames(plotEmat), "LR_reward_raw_zscored", "LR")


# Save as .mat files
writeMat('C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Data/EEG/Export/plotesP3brew_PE_LR.mat',EmatP3brew=plotEmat)


#### late CNV ####

# Topography
data_topo  =  cbind(data.raw, CNV10001500) # the size needs to be exactly the same as the raw eeg data

# delete RTs below 200ms
data_topo$RT = ifelse(data_topo$RT<200,NA,data_topo$RT)

# Rescale the data_topo for the regressions
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_prev[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
# Current trial
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}
for (s in 1:length(unique(data_topo$SubID))){
  data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Feedback locked
for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_efficacy_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(data$SubID))){
  data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]] = scale(data_topo$mbased_reward_feedback_locked[data_topo$SubID == unique(data_topo$SubID)[s]], center = TRUE, scale = TRUE)
}

# Divide by the range to bring to -0.5,0.5 range
data_topo$mbased_efficacy = data_topo$mbased_efficacy/5.837704
data_topo$mbased_efficacy_prev = data_topo$mbased_efficacy_prev/5.837704
data_topo$mbased_efficacy_feedback_locked = data_topo$mbased_efficacy_feedback_locked/5.837704
data_topo$mbased_reward = data_topo$mbased_reward/7.99977
data_topo$mbased_reward_prev = data_topo$mbased_reward_prev/7.99977
data_topo$mbased_reward_feedback_locked = data_topo$mbased_reward_feedback_locked/7.99977

data_topo$feedbackOrder = ifelse(data_topo$feedbackOrder==1,"Efficacy_First","Efficacy_Second") # 1 is first efficacy
data_topo$feedbackOrder = as.factor(data_topo$feedbackOrder)
contrasts(data_topo$feedbackOrder) <- contr.sdif(2)
data_topo$EffLvl = as.factor(data_topo$EffLvl)
data_topo$IsRewarded = as.factor(data_topo$IsRewarded)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
colnames(attr(data_topo$EffLvl, "contrasts")) =  c("Eff_min_NEff") # Change the name of the contrasts
contrasts(data_topo$IsRewarded) <- contr.sdif(2)
colnames(attr(data_topo$IsRewarded, "contrasts")) =  c("Rew_min_NRew") # Change the name of the contrasts
data_topo$Congruency = ordered(data_topo$Congruency, levels = c("Congruent", "Neutral", "Incongruent"))  # first contrast is facilitaion and the second one is interference
contrasts(data_topo$Congruency) <- contr.sdif(3)
colnames(attr(data_topo$Congruency, "contrasts")) =  c("Facilitation", "Interference") # Change the name of the contrasts
data_topo$Accuracy = data_topo$Acc
data_topo$Accuracy = ifelse(data_topo$Accuracy==1,"Correct","Incorrect") 
data_topo$Accuracy = as.factor(data_topo$Accuracy)
contrasts(data_topo$Accuracy) <- contr.sdif(2)
colnames(attr(data_topo$Accuracy, "contrasts")) =  c("Inc_min_Corr") # Change the name of the contrasts
data_topo$Rew_eff_fb = as.factor(data_topo$Rew_eff_fb)
data_topo$Rew_eff_fb = ordered(data_topo$Rew_eff_fb, levels = c("NoReward", "Reward", "RewardUnknown"))
contrasts(data_topo$Rew_eff_fb) <- contr.sdif(3)
colnames(attr(data_topo$Rew_eff_fb, "contrasts")) =  c("Rew_min_NRew","RewUnk_min_Rew") # Change the name of the contrasts
data_topo$Eff_rew_fb = as.factor(data_topo$Eff_rew_fb)
data_topo$Eff_rew_fb = ordered(data_topo$Eff_rew_fb, levels = c("NoEfficacy", "Efficacy", "EfficacyUnknown"))
contrasts(data_topo$Eff_rew_fb) <- contr.sdif(3)
colnames(attr(data_topo$Eff_rew_fb, "contrasts")) =  c("Eff_min_NEff","EffUnk_min_Eff") # Change the name of the contrasts
data_topo$LR = scale(data_topo$LR, scale= FALSE, center = TRUE)
contrasts(data_topo$EffLvl) <- contr.sdif(2)
contrasts(data_topo$IsRewarded) <- contr.sdif(2)


# estimates
plotEmat = as.data.frame(matrix(nrow = 65, ncol=3))  #as many rows as the number of the electrodes, ncol is the number of fixed effects in the model

# Set the priors (Based on Fromero et al., 2020)
prior = c(
  prior(normal(-0.16, 1.59), class = Intercept),
  prior(normal(-0.30, 0.73), class = b, coef=mbased_efficacy_prev),
  prior(normal(0, 0.73), class = b, coef=mbased_reward_prev))

# fit a model to each electrode and extract the estimates 
j=1
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependet variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod <- 
    brm(DV ~ mbased_efficacy_prev + mbased_reward_prev + (1|SubID),
        data=data_topo[data_topo$RT>200,],
        family=gaussian(),
        prior = prior,
        iter  = 1000,
        warmup=3000,
        inits = 0)
  
  k <- summary(Plotmod)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
}

# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")

# Save as .mat files
writeMat('C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Data/EEG/Export/ploteslateCNV.mat',EmatlateCNV=plotEmat)
