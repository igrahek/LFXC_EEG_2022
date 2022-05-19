#### About the code ####

# Analysis code for the online interval version of the LFXC task. 
# Code written by: Ivan Grahek & Mahalia Prater Fahey (2020 & 2021)

#### Clear environment and import data ####

# clear the environment
rm(list=ls()) 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms, broom, tidyverse, brmstools, BEST, knitr, here, zoo, magrittr, pracma,xtable, Hmisc, ppcor, lme4,MuMIn,MASS, sjPlot, jtools, lmerTest, sjstats,coefplot,R.matlab,RColorBrewer,cowplot,bayesplot,rstan,RColorBrewer,rstudioapi,corrplot,PCAmixdata,psych,rmarkdown,knitr,kable,mclust)

# set seed
set.seed(42) 

## Set current directory and paths, and other script parameters 

# Paths for the data
DataPath=paste0('../../../../../Data/Experiment2/')
date1='2021-06-02'
date2='2021-07-15'
TaskVersion1="TSS_LFXC20210602"
TaskVersion2="TSS_LFXC20210714"


# Criteria for Exclusion 
thresh_minIntervals_performed = 288  # Must have performed the whole experiment 
thresh_maxRT_performed = 2000       # Max RT (ms)
thresh_minRT_performed = 250        # Min RT (ms)
thresh_meanAcc_subj = .1
thresh_quiz_mistakes = 10           # Exclude if there are more than 10 errors for any of the quizzes
complete_Task= 1                  #Exclude if they did not complete the task
thresh_maxTimeonTask = 7800 # Exclude if they took too long to complete the task (in seconds - over 130 minutes)

# Data scrubbing variables
task_ISI = 250                      # constant variable set in task

# Set colors for plotting
# To visualize all of the palettes, used: display.brewer.all()
p.colors.hex<-brewer.pal(n = 8, name = 'Paired')


## Import Raw Participant Data

# First dataset
data.prolific1 = read_csv(paste0(DataPath,"prolific_export_TSS_LFXC20210602.csv"))

#Second dataset
data.prolific2 = read_csv(paste0(DataPath,"prolific_export_TSS_LFXC20210714.csv"))

# Merge the datasets
data.prolific = rbind(data.prolific1,data.prolific2)


## Filter Prolific Participant logs by Approved 
data.prolific.filter<- data.prolific %>%
  dplyr::select(participant_id,status)%>%
  dplyr::filter(status=="APPROVED")%>%
  dplyr::select(participant_id)%>%
  distinct(participant_id)%>%
  dplyr::rename(SubID=participant_id)

## Load in datasets

# Load in the behavioral data
# First dataset
df.trial.raw1<-read.csv(paste0(DataPath,'trialdata_TSS_LFXCGarden_',date1,'.csv'))
#Second dataset
df.trial.raw2<-read.csv(paste0(DataPath,'trialdata_TSS_LFXCGarden_',date2,'.csv'))
# Separate uniqueid from version
df.trial.raw2 <- df.trial.raw2 %>%
  separate(col = uniqueid, sep = ":", c("uniqueid","Version"), remove = TRUE) 

#Merge the datasets
df.trial.raw<-rbind(df.trial.raw1,df.trial.raw2)

#Load in quiz data (how many incorrect responses per quiz question). The quiz data are attention checks (Stroop, efficacy manipulations, etc.)
df.quiz1<-read.csv(paste0(DataPath,'questiondata_TSS_LFXCGarden_',date1,'.csv'), header=F,
                   col.names = c('uniqueid','attribute','value','')) %>%
  spread(attribute,value) %>%
  separate(col = uniqueid, sep = ":", c("SubID","Version"), remove = FALSE) %>%
  filter(SubID %in% df.trial.raw$uniqueid) %>% 
  dplyr::select(SubID,incorrect1,incorrect2,incorrect3,incorrect4,incorrect5,EndTask) %>%
  mutate_at(vars(incorrect1:EndTask),as.character,
            vars(incorrect1:EndTask),as.numeric) %>%  #Make sure that these lines actually work, on some R versions/.packages the order of this line abnd the one above is important
  filter(incorrect1 < thresh_quiz_mistakes, incorrect2 < thresh_quiz_mistakes,
         incorrect3 < thresh_quiz_mistakes, incorrect4 < thresh_quiz_mistakes,
         incorrect5 < thresh_quiz_mistakes, EndTask == complete_Task)

df.quiz2<-read.csv(paste0(DataPath,'questiondata_TSS_LFXCGarden_',date2,'.csv'),header=F,
                   col.names = c('uniqueid','attribute','value')) %>%
  spread(attribute,value) %>%
  separate(col = uniqueid, sep = ":", c("SubID","Version"), remove = FALSE) %>%
  filter(SubID %in% df.trial.raw$uniqueid) %>% 
  dplyr::select(SubID,incorrect1,incorrect2,incorrect3,incorrect4,incorrect5,EndTask) %>%
  mutate_at(vars(incorrect1:EndTask),as.character,
            vars(incorrect1:EndTask),as.numeric) %>% #Make sure that these lines actually work, on some R versions/.packages the order of this line abnd the one above is important
  filter(incorrect1 < thresh_quiz_mistakes, incorrect2 < thresh_quiz_mistakes,
         incorrect3 < thresh_quiz_mistakes, incorrect4 < thresh_quiz_mistakes,
         incorrect5 < thresh_quiz_mistakes, EndTask == complete_Task)

#Merge the quiz data
df.quiz<-rbind(df.quiz1,df.quiz2)


#### Clean the trial-level data ####

# Subjects to exclude
sub.trials<-data.frame(table(df.trial.raw$uniqueid))
# Subjects to exclude if they spent more than 2h on the task
sub.exclude.timeOnTask <-data.prolific$participant_id[which(data.prolific$time_taken>thresh_maxTimeonTask)]  

# Copy the probe responses to all trials in an interval
df.trial = df.trial.raw
df.trial$prompt = as.numeric(df.trial.raw$prompt)
df.trial$rating = as.numeric(df.trial.raw$rating)

# save the original number of participants
participantsPre = length(unique(df.trial$uniqueid))

#Format Trial-level data
df.trial <- df.trial %>%
  filter(Version==TaskVersion1|Version==TaskVersion2, uniqueid %in% df.quiz$SubID) %>% #they did one of our two versions and have quiz data
  filter(!uniqueid %in% sub.exclude.timeOnTask)   # Exclude the subjects who took too long to complete the task

# save the number of excluded based on the criteria above
excluded.time.quizzes = participantsPre - length(unique(df.trial$uniqueid)) 

# Rename subjects to normal numbers
df.trial = df.trial %>% group_by(uniqueid) %>% mutate(SubID=cur_group_id()) 

# Add the probe responses
df.trial <- df.trial %>% 
  
  dplyr::group_by(uniqueid,blockNum,intervalNum,intervalType) %>% 
  dplyr::mutate(reward = max(unique(intervalScorePoints),na.rm=TRUE), #this is the reward feedback
                probeType = max(unique(prompt),na.rm=TRUE), #whether this is efficacy (2), reward (1), or no probe (0)
                probeResponse = max(rating,na.rm=TRUE)) %>% ungroup()

# Turn the no probe response intervals into NA      
df.trial$probeResponse <- ifelse(df.trial$probeResponse==-Inf,NA,df.trial$probeResponse)

# Separate the reward probe response
df.trial$RewardProbeResp <- (ifelse(df.trial$probeType==1,df.trial$probeResponse,NA))/100

# Separate the efficacy probe response
df.trial$EfficacyProbeResp <- (ifelse(df.trial$probeType==2,df.trial$probeResponse,NA))/100

# Create a vector of trial number for each subject's session. This should be sequential and provide a trial number across the whole session.  
tmp.trialdata<-df.trial %>% group_by(uniqueid) %>% 
  dplyr::summarise(numtrials=n(), .groups="keep")
for (t in 1:dim(tmp.trialdata)[1]) { 
  Trials = seq(1,tmp.trialdata$numtrials[t])
  if (t==1) {
    TrialSessionNum <- Trials
  }  else {
    TrialSessionNum <- append(TrialSessionNum,Trials)
  }
}

# # Create vector of trial number for each subject's session
# tmp.trialdata<-df.trial %>% group_by(SubID) %>% 
#   dplyr::summarise(numtrials=n(), .groups="keep")
# df.trialNew<-data.frame()
# for (subject in 1:length(unique(tmp.trialdata$SubID))) {
#   #Plotting
#   plot_data = subset(df.trial,df.trial$SubID == unique(tmp.trialdata$SubID)[subject])
#   
#   plot_data<-plot_data  %>% 
#     dplyr::mutate(wholesessionTrialNum = row_number())
#   df.trialNew<-rbind(df.trialNew,plot_data)
#   
# }

# Create additional variables and reorganize experimental factors
df.trial <- df.trial %>% 
  mutate(AccRT=ifelse(hit==1,rt,NaN),
         type=droplevels(factor(type)),
         congruence=as.numeric(type=='congruent'), # 1=congruent, 0=incongruent
         EfficacyLevel=intervalType,
         reward = reward,
         IntervalsPerBlock = max(unique(intervalNum)),
         IntervalSessionNum = IntervalsPerBlock*(blockNum-1) + intervalNum,
         TrialSessionNum = TrialSessionNum,
         scaledIntervalSessionNum = scale(IntervalSessionNum, center = TRUE,scale = TRUE),
         scaledTrialSessionNum = scale(TrialSessionNum, center = TRUE, scale = TRUE),
         scaledIntervalBlockNum = scale(intervalNum, center = TRUE,scale = TRUE),
         scaledTrialIntervalNum = scale(trialNum, center = TRUE, scale = TRUE),
         scaledIntervalLength = scale(intervalLength, center = TRUE, scale = TRUE),
         log10_AccRT=log10(AccRT),
         scaled_log10_AccRT=scale(log10_AccRT, center = TRUE, scale = TRUE),
         lowerWord=tolower(word),
         ErrorType=case_when(
           hit==0 & response!=lowerWord~ "random",
           hit==0 & response==lowerWord ~ "automatic",
           hit==1~ "NoError"))




#### Clean the data frame for interval-level data ####

df.interval <- df.trial %>% 
  group_by(SubID, Version,blockNum,intervalNum,intervalLength,EfficacyLevel) %>% 
  dplyr::summarise(trialsPerInterval = n(),
                   IntervalISI = (trialsPerInterval-1)*task_ISI, # Calculate Total ISI per interval 
                   Interval_Acc=mean(hit,na.rm=T),
                   Interval_sum_Acc=sum(hit,na.rm=T), 
                   Interval_AccRT=mean(AccRT,na.rm=T),
                   Interval_log10_AccRT=mean(log10_AccRT,na.rm=T),
                   Interval_RT=mean(rt,na.rm=T),
                   Interval_Congruence=mean(congruence,na.rm=T), 
                   scaledIntervalSessionNum = mean(scaledIntervalSessionNum),
                   scaledIntervalBlockNum = mean(scaledIntervalBlockNum),
                   scaledIntervalLength = mean(scaledIntervalLength), 
                   .groups="keep",
                   reward = max(as.numeric(reward),na.rm = T),
                   EfficacyProbeResp = ifelse( !all(is.na(EfficacyProbeResp)), max(EfficacyProbeResp, na.rm=T), NA),
                   RewardProbeResp = ifelse( !all(is.na(RewardProbeResp)), max(RewardProbeResp, na.rm=T), NA))  %>% ungroup() %>% 
  
  mutate(Interval_norm_sum_Acc=((Interval_sum_Acc)/(intervalLength-IntervalISI))*1000, 
         scaledIntervalCong = scale(Interval_Congruence, center = TRUE, scale = TRUE), 
         scale_Interval_log10_AccRT=scale(Interval_log10_AccRT, center = TRUE, scale = TRUE)) 

# Calculate the linear extrapolation of the subjective estimates
df.interval <- df.interval %>%
  group_by(SubID) %>%
  mutate(EfficacyProbeRespLin = na.approx(EfficacyProbeResp, na.rm=F),
         RewardProbeRespLin = na.approx(RewardProbeResp, na.rm=F))  %>% ungroup()

# Create the previous efficacy/reward subjective estimates variables
df.interval = ddply(df.interval,.(SubID),transform,
                    RewardProbeRespLin_prev = append(RewardProbeRespLin,NA,after=0)[-(length(RewardProbeRespLin)+1)],   
                    EfficacyProbeRespLin_prev = append(EfficacyProbeRespLin,NA,after=0)[-(length(EfficacyProbeRespLin)+1)])


# Add the overall Interval number
df.interval = ddply(df.interval,.(SubID),plyr::mutate,
                    Interval = 1:(length(df.interval$intervalNum)/length(unique(df.interval$SubID))))

# Recode efficacy level to numeric for calculating the running average
df.interval = df.interval %>%
  mutate_at("EfficacyLevel",funs(recode(.,"random_low" = "0","performance_low" = "1")))

# Turn into numeric
df.interval$EfficacyLevel = as.numeric(as.character(df.interval$EfficacyLevel))
df.interval$reward = as.numeric(as.character(df.interval$reward))

# Range normalize reward magnitude to 0-1
df.interval <- df.interval %>% 
  group_by(SubID) %>% 
  mutate(reward=(reward-min(reward,na.rm = T))/(max(reward,na.rm = T)-min(reward,na.rm = T)))%>% ungroup()


# Create previous efficacy variable
df.interval = ddply(df.interval,.(SubID),transform,
                    Efficacy_T0 = EfficacyLevel,   
                    Efficacy_T1 = append(EfficacyLevel,NA,after=0)[-(length(EfficacyLevel)+1)],
                    Efficacy_T2 = append(EfficacyLevel,c(NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2))],
                    Efficacy_T3 = append(EfficacyLevel,c(NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3))],
                    Efficacy_T4 = append(EfficacyLevel,c(NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4))],
                    Efficacy_T5 = append(EfficacyLevel,c(NA,NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4),(length(EfficacyLevel)+5))],
                    Efficacy_T6 = append(EfficacyLevel,c(NA,NA,NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4),(length(EfficacyLevel)+5),(length(EfficacyLevel)+6))],
                    Efficacy_T7 = append(EfficacyLevel,c(NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4),(length(EfficacyLevel)+5),(length(EfficacyLevel)+6),(length(EfficacyLevel)+7))],
                    Efficacy_T8 = append(EfficacyLevel,c(NA,NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4),(length(EfficacyLevel)+5),(length(EfficacyLevel)+6),(length(EfficacyLevel)+7),(length(EfficacyLevel)+8))],
                    Efficacy_T9 = append(EfficacyLevel,c(NA,NA,NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(EfficacyLevel)+1),(length(EfficacyLevel)+2),(length(EfficacyLevel)+3),(length(EfficacyLevel)+4),(length(EfficacyLevel)+5),(length(EfficacyLevel)+6),(length(EfficacyLevel)+7),(length(EfficacyLevel)+8),(length(EfficacyLevel)+9))])


# Create previous reward variable
df.interval = ddply(df.interval,.(SubID),transform,
                    Reward_T0 = reward,   
                    Reward_T1 = append(reward,NA,after=0)[-(length(reward)+1)],
                    Reward_T2 = append(reward,c(NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2))],
                    Reward_T3 = append(reward,c(NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3))],
                    Reward_T4 = append(reward,c(NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4))],
                    Reward_T5 = append(reward,c(NA,NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4),(length(reward)+5))],
                    Reward_T6 = append(reward,c(NA,NA,NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4),(length(reward)+5),(length(reward)+6))],
                    Reward_T7 = append(reward,c(NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4),(length(reward)+5),(length(reward)+6),(length(reward)+7))],
                    Reward_T8 = append(reward,c(NA,NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4),(length(reward)+5),(length(reward)+6),(length(reward)+7),(length(reward)+8))],
                    Reward_T9 = append(reward,c(NA,NA,NA,NA,NA,NA,NA,NA,NA),after=0)[-c((length(reward)+1),(length(reward)+2),(length(reward)+3),(length(reward)+4),(length(reward)+5),(length(reward)+6),(length(reward)+7),(length(reward)+8),(length(reward)+9))])


# Turn categorical variables into factors
# Recode and turn into a factor
EfficacyLag = c("Efficacy_T0",
                "Efficacy_T1",
                "Efficacy_T2",
                "Efficacy_T3",
                "Efficacy_T4",
                "Efficacy_T5",
                "Efficacy_T6",
                "Efficacy_T7",
                "Efficacy_T8",
                "Efficacy_T9")

# Recode and turn into a factor
RewardLag = c("Reward_T0",
              "Reward_T1",
              "Reward_T2",
              "Reward_T3",
              "Reward_T4",
              "Reward_T5",
              "Reward_T6",
              "Reward_T7",
              "Reward_T8",
              "Reward_T9")
RunAvgWindow=5

#Add the running averages for reward, efficacy, and performance
df.interval =  df.interval %>%
  dplyr::group_by(SubID)%>%
  dplyr::mutate(
    runAvgEfficacy = rollapply(Efficacy_T0,list(seq(-RunAvgWindow, -1)),mean, partial=TRUE,fill = NA, na.rm = TRUE, align = "right"),
    runAvgReward = rollapply(Reward_T0,list(seq(-RunAvgWindow, -1)),mean, partial=TRUE,fill = NA, na.rm = TRUE, align = "right"),
    runAvgPerformance = rollapply(Interval_norm_sum_Acc,list(seq(-3, -1)),mean, partial=TRUE,fill = NA, na.rm = TRUE, align = "right"))




#### Import the learning rates for efficacy ###### 


learning_rates = read.csv("../../../RL_fitting/4_Intercept_Pos_and_neg_learning_rate/results/learning_rates_efficacy.csv", header = F)

# add the subject names
learning_rates$SubID = unique(df.interval$SubID)

# rename the variable names
colnames(learning_rates)[1] <- "positive_learning_rate"
colnames(learning_rates)[2] <- "negative_learning_rate"
colnames(learning_rates)[3] <- "initial_bias"
colnames(learning_rates)[5] <- "subject"
colnames(learning_rates)[6] <- "BIC"

learning_rates_efficacy = learning_rates

# add a new variable in the main dataset in which efficacy is coded as 1 or 0
df.interval$EffLvl_forLR = df.interval$EfficacyLevel

#### Calculate the model-based efficacy estimate ######

# for each subject
for (s in 1:length(unique(df.interval$SubID))) {
  # for each trial
  for (t in 1:length(unique(df.interval$Interval))) {
    # if this is the first trial use the inital estimate
    if (t == 1) {
      v = learning_rates$initial_bias[s]
    } else {
      v = v
    }
    
    # save the BIC value
    df.interval$BIC_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                               df.interval$Interval == t] = learning_rates$BIC[s]
    # calculate the signed prediction error for this subject for this trial (difference between the actual and the expected efficacy)
    df.interval$signed_PE_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                     df.interval$Interval == t] = df.interval$EffLvl_forLR[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                             df.interval$Interval == t] - v
    
    # calculate the unsigned prediction error for this subject for this trial (absolute difference between the actual and the expected efficacy)
    df.interval$unsigned_PE_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                       df.interval$Interval == t] = abs(df.interval$EffLvl_forLR[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                   df.interval$Interval == t] - v)
    
    # save the difference between the positive and the negative LR
    df.interval$LR_efficacy_Pos_min_Neg[df.interval$SubID == unique(df.interval$SubID)[s] &
                                          df.interval$Interval == t] = learning_rates$positive_learning_rate[s] - learning_rates$negative_learning_rate[s]
    
    # if the signed pe is positive
    if (df.interval$signed_PE_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                       df.interval$Interval == t] > 0) {
      
      
      
      # calculate the size of the update (the learning rate times the prediction error)
      df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                           df.interval$Interval == t] = learning_rates$positive_learning_rate[s] * df.interval$signed_PE_efficacy[df.interval$SubID ==
                                                                                                                                                    unique(df.interval$SubID)[s] & df.interval$Interval == t]
      # save the absolute value of the update
      df.interval$mbased_efficacy_absolute_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                    df.interval$Interval == t] =  abs(df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                                           df.interval$Interval == t])
      
      # calculate the effiacy estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      df.interval$mbased_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                    df.interval$Interval == t] = v + df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # update the value
      v = v + df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # save the learning rate value
      df.interval$LR_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                df.interval$Interval == t] = learning_rates$positive_learning_rate[s]
      
      # save the positive learning rate value
      df.interval$LR_efficacy_positive[df.interval$SubID == unique(df.interval$SubID)[s] &
                                         df.interval$Interval == t] = learning_rates$positive_learning_rate[s]
      
      # if the rpe is negative
    } else {
      
      # calculate the size of the update (the learning rate times the prediction error)
      df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                           df.interval$Interval == t] = learning_rates$negative_learning_rate[s] * df.interval$signed_PE_efficacy[df.interval$SubID ==
                                                                                                                                                    unique(df.interval$SubID)[s] & df.interval$Interval == t]
      # save the absolute value of the update
      df.interval$mbased_efficacy_absolute_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                    df.interval$Interval == t] =  abs(df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                                           df.interval$Interval == t])
      
      # calculate the effiacy estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      df.interval$mbased_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                    df.interval$Interval == t] = v + df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # update the value
      v = v + df.interval$mbased_efficacy_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # save the learning rate value
      df.interval$LR_efficacy[df.interval$SubID == unique(df.interval$SubID)[s] &
                                df.interval$Interval == t] = learning_rates$negative_learning_rate[s]
      
      # save the positive learning rate value
      df.interval$LR_efficacy_negative[df.interval$SubID == unique(df.interval$SubID)[s] &
                                         df.interval$Interval == t] = learning_rates$negative_learning_rate[s]
      
    }
    
    
  }
}


#### Import the learning rates for reward rate###### 
learning_rates = read.csv("../../../RL_fitting/4_Intercept_Pos_and_neg_learning_rate/results/learning_rates_reward.csv", header = F)

# add the subject names
learning_rates$SubID = unique(df.interval$SubID)

# rename the variable names
colnames(learning_rates)[1] <- "positive_learning_rate"
colnames(learning_rates)[2] <- "negative_learning_rate"
colnames(learning_rates)[3] <- "initial_bias"
colnames(learning_rates)[5] <- "subject"
colnames(learning_rates)[6] <- "BIC"

learning_rates_reward = learning_rates


# add a new variable in the main dataset in which efficacy is coded as 1 or 0
df.interval$IsRewarded_forLR = df.interval$reward

#### Calculate the model-based reward rate estimate ######

# for each subject
for (s in 1:length(unique(df.interval$SubID))) {
  # for each trial
  for (t in 1:length(unique(df.interval$Interval))) {
    # if this is the first trial use the inital estimate
    if (t == 1) {
      v = learning_rates$initial_bias[s]
    } else {
      v = v
    }
    
    # save the BIC value
    df.interval$BIC_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                             df.interval$Interval == t] = learning_rates$BIC[s]
    
    # calculate the prediction error for this subject for this trial (difference between the actual and the expected efficacy)
    df.interval$signed_PE_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                                   df.interval$Interval == t] = df.interval$IsRewarded_forLR[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                               df.interval$Interval == t] - v
    
    # calculate the unsigned prediction error for this subject for this trial (absolute difference between the actual and the expected efficacy)
    df.interval$unsigned_PE_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                                     df.interval$Interval == t] = abs(df.interval$IsRewarded_forLR[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                     df.interval$Interval == t] - v)
    
    # if the rpe is positive
    if (df.interval$signed_PE_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                                     df.interval$Interval == t] > 0) {
      
      # calculate the size of the update (the learning rate times the prediction error)
      df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                         df.interval$Interval == t] = learning_rates$positive_learning_rate[s] * df.interval$signed_PE_reward[df.interval$SubID ==
                                                                                                                                                unique(df.interval$SubID)[s] & df.interval$Interval == t]
      # save the absolute value of the update
      df.interval$mbased_reward_absolute_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                  df.interval$Interval == t] =  abs(df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                                       df.interval$Interval == t])
      
      # calculate the reward estimate (expected reward plus the update - the learning rate-weighted prediction error)
      df.interval$mbased_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                                  df.interval$Interval == t] = v + df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # update the value
      v = v + df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # save the learning rate value
      df.interval$LR_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                              df.interval$Interval == t] = learning_rates$positive_learning_rate[s]
      
      # save the positive learning rate value
      df.interval$LR_reward_positive[df.interval$SubID == unique(df.interval$SubID)[s] &
                                       df.interval$Interval == t] = learning_rates$positive_learning_rate[s]
      
      # if the rpe is negative
    } else {
      
      # calculate the size of the update (the learning rate times the prediction error)
      df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                         df.interval$Interval == t] = learning_rates$negative_learning_rate[s] * df.interval$signed_PE_reward[df.interval$SubID ==
                                                                                                                                                unique(df.interval$SubID)[s] & df.interval$Interval == t]
      # save the absolute value of the update
      df.interval$mbased_reward_absolute_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                  df.interval$Interval == t] =  abs(df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] &
                                                                                                                       df.interval$Interval == t])
      
      # calculate the reward estimate (expected efficacy plus the update - the learning rate-weighted prediction error)
      df.interval$mbased_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                                  df.interval$Interval == t] = v + df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # update the value
      v = v + df.interval$mbased_reward_update[df.interval$SubID == unique(df.interval$SubID)[s] & df.interval$Interval == t]
      
      # save the learning rate value
      df.interval$LR_reward[df.interval$SubID == unique(df.interval$SubID)[s] &
                              df.interval$Interval == t] = learning_rates$negative_learning_rate[s]
      
      # save the positive learning rate value
      df.interval$LR_reward_negative[df.interval$SubID == unique(df.interval$SubID)[s] &
                                       df.interval$Interval == t] = learning_rates$negative_learning_rate[s]
      
    }
    
    
  }
}


# Create a previous value variable (for predicting everything before the new feedback participants rely on the value from the previous trial)
df.interval = ddply(df.interval,.(SubID),transform,
                    mbased_reward_prev = append(mbased_reward,NA,after=0)[-(length(mbased_reward)+1)],
                    mbased_efficacy_prev = append(mbased_efficacy,NA,after=0)[-(length(mbased_efficacy)+1)])


#### Exclude too fast and too slow RTs (<250ms)####

# Too fast
df.trial <- df.trial %>%
  filter(rt>thresh_minRT_performed)

# Too slow
df.trial <- df.trial %>%
  filter(rt<thresh_maxRT_performed)


#### Exclude based on accuracy ####

summary_table = ddply(df.trial,.(SubID),plyr::summarize,
                      MeanRT=mean(AccRT,na.rm=TRUE), # mean RT per condition
                      SDRT=sd(AccRT,na.rm=TRUE),
                      MeanAcc=mean(hit,na.rm=TRUE))

print(kable(summary_table, align="l", caption="Accuracy per subject"))

sub.exclude.accuracy = summary_table$SubID[which(summary_table$MeanAcc<0.7)]

# #Exclude from interval level data
df.interval <- df.interval %>%
  filter(!SubID %in% sub.exclude.accuracy)

#Exclude from trial level data
df.trial <- df.trial %>%
  filter(!SubID %in% sub.exclude.accuracy)


#### Exclude participants based on efficacy learning ####


## Clustering

# Exclude subjects with flat efficacy estimates
LR_all = cbind(learning_rates_efficacy,learning_rates_reward)

colnames(LR_all) = c("PosLR_Eff",
                     "NegLR_Eff",
                     "Bias_Eff",
                     "Variance_Eff",
                     "Subject",
                     "BIC_Eff",
                     "SubID",
                     "PosLR_Rew",
                     "NegLR_Rew",
                     "Bias_Rew",
                     "Variance_Rew",
                     "Subject",
                     "BIC_Rew",
                     "SubID2")
LR_all = subset(LR_all,select=-c(Subject,SubID,SubID2))

# Exclude the people with low accuracy
LR_all <- LR_all %>%
  filter(!Subject %in% sub.exclude.accuracy)

# Take only the Learning rates data for efficacy
LR = subset(LR_all,select=c(PosLR_Eff,NegLR_Eff))
clPairs(LR)

# BICs for the models
BIC <- mclustBIC(LR)
# plot(BIC)
# summary(BIC)

# Fit the models
mod1 <- Mclust(LR, x = BIC)
summary(mod1, parameters = TRUE)

#Plot the model
pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure_S5_B.pdf")
p=plot(mod1, what = "classification")

p

# Add Classification to the data
LR_all$Class_Eff = mod1[["classification"]]

x = (subset(LR_all,LR_all$Class_Eff==1))
length(x$PosLR_Eff)

# Plot the positive LR histograms
LR_all$Class_Eff=as.factor(LR_all$Class_Eff)

pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS5_C1.pdf")

p=ggplot(LR_all, aes(x=PosLR_Eff,fill=Class_Eff,color=Class_Eff)) + 
  geom_histogram( position="identity", alpha=0.5,bins=80)+
  ggtitle('Positive learning rates') + 
  theme_classic(base_size = 15) 
p
dev.off()
p

pdf(file="C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS5_C2.pdf")

p=ggplot(LR_all, aes(x=NegLR_Eff,fill=Class_Eff,color=Class_Eff)) + 
  geom_histogram( position="identity", alpha=0.5,bins=80)+
  ggtitle('Negative learning rates') + 
  theme_classic(base_size = 15) 

p
dev.off()
p

# Add Classification to the data
LR_all$Class_Eff = mod1[["classification"]]

sub.exclude.efficacy = subset(LR_all,LR_all$Class_Eff==1)
sub.exclude.efficacy = sub.exclude.efficacy$Subject

# Add one subject who has 0 variance in the reward probe responses, thus having NAs after the within-subject z-scoring
sub.exclude.efficacy = append(sub.exclude.efficacy,40)

#Exclude from interval level data
df.interval <- df.interval %>%
  filter(!SubID %in% sub.exclude.efficacy)  

#Exclude from trial level data
df.trial <- df.trial %>%
  filter(!SubID %in% sub.exclude.efficacy) 



#### Join the efficacy estimates with the trial-data ####

df.trial.combined<-df.trial

df.trial.combined = join(df.trial.combined,df.interval, by = c("SubID", "blockNum", "intervalNum"),type="left")

# Change back the names of the columns
names(df.trial.combined)[names(df.trial.combined) == "scaledIntervalLength.x"] <- "scaledIntervalLength"
names(df.trial.combined)[names(df.trial.combined) == "scaledIntervalSessionNum.x"] <- "scaledIntervalSessionNum"

df.interval.combined<-df.interval


#### Re-scale the data ####


# Center the continuous variables 
# df.interval.combined$EfficacyProbeResp = scale(df.interval.combined$EfficacyProbeResp, scale= FALSE, center = TRUE)
df.interval.combined$EfficacyProbeRespLin = scale(df.interval.combined$EfficacyProbeRespLin, scale= FALSE, center = TRUE)
df.interval.combined$EfficacyProbeRespLin_prev = scale(df.interval.combined$EfficacyProbeRespLin_prev, scale= FALSE, center = TRUE)

# Z-score within subject
for (s in 1:length(unique(df.interval.combined$SubID))){
  df.interval.combined$mbased_efficacy_prev[df.interval.combined$SubID == unique(df.interval.combined$SubID)[s]] = scale(df.interval.combined$mbased_efficacy_prev[df.interval.combined$SubID == unique(df.interval.combined$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(df.interval.combined$SubID))){
  df.interval.combined$mbased_reward_prev[df.interval.combined$SubID == unique(df.interval.combined$SubID)[s]] = scale(df.interval.combined$mbased_reward_prev[df.interval.combined$SubID == unique(df.interval.combined$SubID)[s]], center = TRUE, scale = TRUE)
}

# Bring to the 0-1 range
# Bring to the 0-1 range
df.interval.combined$mbased_efficacy_prev = df.interval.combined$mbased_efficacy_prev/(range(df.interval.combined$mbased_efficacy_prev,na.rm = T)[2]-range(df.interval.combined$mbased_efficacy_prev,na.rm = T)[1])
df.interval.combined$mbased_reward_prev = df.interval.combined$mbased_reward_prev/(range(df.interval.combined$mbased_reward_prev,na.rm = T)[2]-range(df.interval.combined$mbased_reward_prev,na.rm = T)[1])


df.interval.combined$Interval = scale(df.interval.combined$Interval,scale= TRUE, center = TRUE)

# Center the continuous variables 
df.trial.combined$EfficacyProbeRespLin = scale(df.trial.combined$EfficacyProbeRespLin, scale= FALSE, center = TRUE)
df.trial.combined$EfficacyProbeRespLin_prev = scale(df.trial.combined$EfficacyProbeRespLin_prev, scale= FALSE, center = TRUE)

# Z-score within subject
for (s in 1:length(unique(df.trial.combined$SubID))){
  df.trial.combined$mbased_efficacy_prev[df.trial.combined$SubID == unique(df.trial.combined$SubID)[s]] = scale(df.trial.combined$mbased_efficacy_prev[df.trial.combined$SubID == unique(df.trial.combined$SubID)[s]], center = TRUE, scale = TRUE)
}

for (s in 1:length(unique(df.trial.combined$SubID))){
  df.trial.combined$mbased_reward_prev[df.trial.combined$SubID == unique(df.trial.combined$SubID)[s]] = scale(df.trial.combined$mbased_reward_prev[df.trial.combined$SubID == unique(df.trial.combined$SubID)[s]], center = TRUE, scale = TRUE)
}
# Bring to the 0-1 range
df.trial.combined$mbased_efficacy_prev = df.trial.combined$mbased_efficacy_prev/(range(df.trial.combined$mbased_efficacy_prev,na.rm = T)[2]-range(df.trial.combined$mbased_efficacy_prev,na.rm = T)[1])
df.trial.combined$mbased_reward_prev = df.trial.combined$mbased_reward_prev/(range(df.trial.combined$mbased_reward_prev,na.rm = T)[2]-range(df.trial.combined$mbased_reward_prev,na.rm = T)[1])


df.trial.combined$Interval = scale(df.trial.combined$Interval,scale= TRUE, center = TRUE)
df.trial.combined$runAvgEfficacy = scale(df.trial.combined$runAvgEfficacy, scale= FALSE, center = TRUE)
df.trial.combined$runAvgReward = scale(df.trial.combined$runAvgReward, scale= FALSE, center = TRUE)
df.trial.combined$runAvgPerformance = scale(df.trial.combined$runAvgPerformance, scale= FALSE, center = TRUE)

# Change the name of the congruency variable
df.trial.combined$Congruency = df.trial.combined$type

# Change accuracy to numeric to fit
df.trial.combined$hit = as.numeric(df.trial.combined$hit)






#### Save out for brms ####
write.csv(df.trial.combined,"Exp2Data_0.1Acc.csv", row.names=FALSE)


#### Pre-processing summary ####
preproSummary = as.data.frame(matrix())
preproSummary$V1 = "NSubs"

# Time on task and quizzes
preproSummary$Initial = length(unique(df.trial.raw$uniqueid))

# Efficacy clustering
preproSummary$EfficacyLearnigRates = length(sub.exclude.efficacy)  

# Accuracy
preproSummary$Accuracy = length(sub.exclude.accuracy) 

# Time on task and quizzes
preproSummary$Quizzes_Time = excluded.time.quizzes 

# Demographics
final_participants <- data.prolific %>%
  filter(participant_id %in% unique(df.trial$uniqueid))

# Median age
preproSummary$MedianAge = median(final_participants$age)

# Gender
preproSummary$Male = length(final_participants[ which(final_participants$Sex=='Male'),]$Sex)
preproSummary$Female = length(final_participants[ which(final_participants$Sex=='Female'),]$Sex)


# Final
preproSummary$Final = length(unique(df.trial.combined$uniqueid)) 

write.csv(preproSummary, file = "PreprocessingSummary_0.1Acc.csv")






