#!/usr/bin/env bash
#### Bradley 2023
# Submit scripts for chromosome 1 only
# bsub -o logs/005_006_master-%J-output.log -e logs/005_006_master-%J-error.log -q normal -G team152 -n 1 -M 5000 -a "memlimit=True" -R "select[mem>5000] rusage[mem=5000] span[hosts=1]" -J "005_006_master" < run_scripts/005_006_submit_chrom1.sh 


# Make sure am in the script dir
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE

# Submit array for chromosome 1 genes (Secretory, TI)
bsub -o logs/saige_array_test-%J-%I-output.log -e logs/saige_array_test-%J-%I-error.log -q normal -G team152 -n 1 -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -J "saige_array_chrom1[1-966]%200" < testing_scripts/005-run_SAIGE_1_2_3_chrom1.sh 

# Submit job to tidy this up and optimise PCs
# bsub -o logs/saige_clean_chrom1-%J-output.log -e logs/saige_clean_chrom1-%J-error.log -q normal -G team152 -n 1 -M 50000 -a "memlimit=True" -R "select[mem>50000] rusage[mem=50000] span[hosts=1]" -J "saige_clean_chrom1" < testing_scripts/006-summarise_SAIGE_1_2_3_chrom1.sh 

