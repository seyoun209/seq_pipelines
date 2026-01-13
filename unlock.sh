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
        '01_run_RNA_preprocess.sbatch')
            snakemake -j 1 --unlock -s workflow/rna_process/01_rna_preprocess.smk --configfile config/RNA_config.yaml --profile config/profile_slurm_rna
            ;;
	'02_RNAmergeSignal.sbatch')
		snakemake -j 1 --unlock -s workflow/rna_process/02_mergeSignal.smk --configfile config/RNA_config.yaml --profile config/profile_slurm_rna
		;;
        '03_RNA_rmats.sbatch')
                snakemake -j 1 --unlock -s workflow/rna_process/03_rmats.smk --configfile config/RNA_config.yaml --profile config/profile_slurm_rna
                ;;
    '01_02_run_atac_peakcall.sbatch')
            snakemake -j 1 --unlock -s workflow/atac_process/02_02_wasp_callingPeak.smk  --configfile config/ATAC_config.yaml --profile config/profile_slurm_atac
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
