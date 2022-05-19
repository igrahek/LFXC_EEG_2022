#### About the code ####

# Fit the RL models in Stan

# =============================================================================
#### Import the data #### 
# =============================================================================

# clear the environment
rm(list=ls()) 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,tidyverse,rstan)

# set seed
set.seed(42) 

## Import the data
data = fread('../data/RL_fit_data_Experiment1.csv')

## Put the data into the Stan format

# The number of subjects
nSubjects = length(unique(data$SubID))

# Add trial number
data$Trial = rep(1:288,nSubjects)

# Add trial number
nTrials = length(unique(data$Trial))

# Rename subjects to normal numbers
data$SubID = rep(1:nSubjects,each=nTrials)

# Replace the NAs with 9999 (Stan doesn't process NAs)
data$EfficacyProbeResp = ifelse(is.na(data$EfficacyProbeResp)==1,9999,data$EfficacyProbeResp)

dataList <- list(nSubjects=nSubjects,
                 nTrials=nTrials, 
                 ProbeResp = t(as.matrix(lapply(split(data,by="SubID"), function(x) x$EfficacyProbeResp) %>% bind_rows())),
                 Feedback = t(as.matrix(lapply(split(data,by="SubID"), function(x) x$EffLvl) %>% bind_rows())))
                 
                 
# =============================================================================
#### Running Stan #### 
# =============================================================================
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

modelFile = '../models/2LR_model.stan'

nIter     <- 10000
nChains   <- 4 
nWarmup   <- 8000
nThin     <- 1

cat("Estimating", modelFile, "model... \n")
startTime = Sys.time(); print(startTime)
cat("Calling", nChains, "simulations in Stan... \n")

fit_rl <- stan(modelFile, 
               data    = dataList, 
               chains  = nChains,
               iter    = nIter,
               warmup  = nWarmup,
               thin    = nThin,
               init    = "random",
               seed    = 1450154626)

cat("Finishing", modelFile, "model simulation ... \n")
endTime = Sys.time(); print(endTime)  
cat("It took",as.character.Date(endTime - startTime), "\n")

saveRDS(fit_rl,file="../output/2LREfficacy.rds")


# Read in the model
# fit_summary=readRDS("2LR_model.rds")

# Summary
fit_summary <- summary(fit_rl)$summary

fit_summary = as.data.table(fit_summary[,c("mean")])


pos_LRs = as.numeric(fit_summary[5:44,])
neg_LRs = as.numeric(fit_summary[45:84,])

# scatter
plot(pos_LRs$V1,neg_LRs$V1)

# histograms
hist(pos_LRs$V1,breaks = 40)
hist(neg_LRs$V1,breaks = 40)
