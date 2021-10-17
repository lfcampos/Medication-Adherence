#!/bin/bash
#SBATCH -n 4 													# Number of cores requested
#SBATCH -N 1 													# Ensure that all cores are on one machine
#SBATCH -t 7200													# Runtime in minutes
#SBATCH -p shared												# Partition to submit to
#SBATCH --mem-per-cpu=1000 										# Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o /n/home01/khunter33/Medication-Adherence/log/%j.out 	# Standard out goes to this file
#SBATCH -e /n/home01/khunter33/Medication-Adherence/log/%j.err 	# Standard err goes to this filehostname
#SBATCH --mail-type=ALL        									# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kristenbhunter@gmail.com 					# Email

source new-modules.sh
module load gcc/8.2.0-fasrc01
module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
R CMD BATCH --quiet --no-restore --no-save /n/home01/khunter33/Medication-Adherence/predict/predict_onestep.R /n/home01/khunter33/Medication-Adherence/log/predict_$(date +"%Y%m%d_%H%M").out



