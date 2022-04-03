#!/bin/sh
 # computing info
#SBATCH --job-name=simu
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ERROR,ARRAY_TASKS
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=7000m 
#SBATCH --time=18:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log

export out_file=/home/dyarger/multivariate_matern/out_simu.out

# The application(s) to execute along with its input arguments and options:
cd /home/dyarger/multivariate_matern/
module load R/4.1.0
R CMD BATCH code/simulation_real.R $out_file
