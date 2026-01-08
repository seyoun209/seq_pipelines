#! /bin/bash -login
#SBATCH -J rna_preprocess
#SBATCH -t 11-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2  
#SBATCH -p general
#SBATCH --mem=4gb
#SBATCH -o workflow_logs/%x-%j.out

## Exit if any command fails
set -e

## Load required modules
module load python/3.12.4

## Create and activate virtual environment with requirements
python3 -m venv env
source env/bin/activate && pip3 install -r config/requirements.txt
#pip install --upgrade pip

## Make directory for slurm logs
mkdir -p rna_output/logs_slurm workflow_logs

#Run
snakemake -s workflow/rna_process/01_rna_preprocess.smk \
  --configfile config/RNA_config.yaml \
  --profile config/profile_slurm_rna \
  -j 120 --max-jobs-per-second 10 --slurm-init-seconds-before-status-checks 30 --max-status-checks-per-second 0.5 \
  --rerun-incomplete -p --latency-wait 150




