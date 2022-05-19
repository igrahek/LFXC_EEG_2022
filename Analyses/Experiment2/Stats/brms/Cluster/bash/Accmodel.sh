#!/bin/bash

#SBATCH --account=carney-ashenhav-condo
#SBATCH --time=400:00:00
#SBATCH --mem=48G
#SBATCH -n 1
#SBATCH -c 5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ivan_grahek@brown.edu
#SBATCH -J Accmodel
#SBATCH -o R-%x.%j.out
module load R/4.1.0 gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018 v8/3.14.5
R CMD BATCH --no-save --quiet --slave ../models/AccModel.R 
