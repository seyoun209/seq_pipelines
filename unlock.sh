#! /bin/bash -login

## Exit if any command fails
set -e

## Load required modules
module load python/3.12.4

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecify which workflow to unlock (i.e. run_RNAprocessing)'
            exit 2
            ;;
        '01_01_run_atac_core.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/00_atac_preprocess.smk --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
    '01_02_run_atac_peakcall.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/02_02_wasp_callingPeak.smk  --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
