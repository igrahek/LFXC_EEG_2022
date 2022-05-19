#### About the code ####

# Experiment: LFXC_EEG 
# Code written by: Ivan Grahek (2018-2021)
# Description: Code for the analysis of behavioral data for the LFXC_EEG project.  

#### Clear the environment and import the data ####---------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstan,plyr,Rmisc,reshape2,brms, broom, tidyverse, knitr, here, zoo, magrittr, pracma,xtable, Hmisc, ppcor,MuMIn,MASS, sjPlot, jtools, lmerTest, sjstats,coefplot,R.matlab,RColorBrewer,ggeffects,glmmTMB, interactions)

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


#### brms settings ####


#help stan run faster
# rstan_options(auto_write = TRUE)
cores=options(mc.cores = parallel::detectCores())
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')


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

# Take only the correct and not too fast responses
data_topo = data_topo[data_topo$Acc == 1 & data_topo$RT>200,]

# Create a list of data
data_topo_list = list()

# Set the working directory where to save the models
setwd(here("Analyses/Stats/brms/brms_models"))


# fit a model to each electrode and extract the estimates 
j=1
for (i in (ncol(data.raw)+1):(ncol(data_topo))){  #go through each electrode and make it the dependent variable
  data_topo$DV  <- data_topo[,i]
  
  # fit a model for each electrode
  Plotmod <- 
    brm(DV ~ EffLvl * mbased_efficacy_prev + mbased_reward_feedback_locked + (1|SubID),
        data=data_topo[data_topo$Acc == 1 & data_topo$RT>200,],
        family=gaussian(),
        prior = prior,
        iter  = 10000,
        warmup=8000)
  
  k <- summary(Plotmod)
  if (i ==(ncol(data.raw)+1)){
    f <- dimnames(k$fixed)[[1]]
    colnames(plotEmat) <-  f}
  plotEmat[j,] <- k$fixed[,1] #estimate
  j=j+1
}

# fit a model to each electrode and extract the estimates 


# Change the stuff in the variable names which Matlab doesn't like
colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "mbased_reward_feedback_locked", "Rew_estimate")
colnames(plotEmat) <- sedit(colnames(plotEmat), "EffLvlEff_min_NEffbymbased_efficacy_prev", "Eff_min_NEffbyEff_estimate")

# Save as .mat files
writeMat('../output/plotesP3beff.mat',EmatP3beff=plotEmat)


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
        iter  = 10000,
        warmup=8000)
  
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
writeMat('../output/plotesP3beff_PE_LR.mat',EmatP3beff=plotEmat)



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
        iter  = 10000,
        warmup=8000)
  
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
writeMat('../output/plotesP3brew.mat',EmatP3brew=plotEmat)


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
        iter  = 10000,
        warmup=8000)
  
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
writeMat('../output/plotesP3brew_PE_LR.mat',EmatP3brew=plotEmat)


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

# Set the priors (Based on Fromero et al., 2021)
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
        iter  = 10000,
        warmup=8000)
  
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
writeMat('../output/ploteslateCNV.mat',EmatlateCNV=plotEmat)
