#!/bin/bash

#SBATCH --account=carney-ashenhav-condo
#SBATCH --time=200:00:00
#SBATCH --mem=48G
#SBATCH -n 1
#SBATCH -c 5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ivan_grahek@brown.edu
#SBATCH -J Experiment2LFXCModel1
#SBATCH -o R-%x.%j.out

python ../../script/model_fitting.py ../../models/LFXC/Model1.json
python ../../script/model_posterior_predictive_check.py LFXC 0 Model1 data_hddm_LFXC_interval_allTrials.csv 500

