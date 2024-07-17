#!/usr/bin/env bash
#### Bradley 2023
#### Generating summary plots on the cis- and trans-eQTL tests. This also includes a comparison against the pseudo-bulk results
# bsub -o logs/saige_summary-%J-output.log -e logs/saige_summary-%J-error.log -q long -G team152 -n 1 -M 120000 -a "memlimit=True" -R "select[mem>120000] rusage[mem=120000] span[hosts=1]" -J "saige_summary" < testing_scripts/016-Summary_plots.sh 

module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
level="T_cell_CD8_1"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="label__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=1
condition_col="NULL"
condition="NULL"
covariates="age_imputed,sex"
covariates_cell=""
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5
use_GRM=FALSE
repo_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE"

# Construct the category directory path
if [ "$condition_col" != "NULL" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi
echo "catdir is ${catdir}"

# Load n_expr_pcs
optim_npcs_file=${catdir}/run_params/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Execute python scripts
python ${repo_dir}/testing_scripts/bin/summary_cis_trans_plots.py \
        --catdir $catdir \
        --chromosomes_cis "1-22" \
        --pb_dir "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123" \
        --phenotype__file $phenotype__file \
        --n_expr_pcs $n_expr_pcs \
        --general_file_dir $general_file_dir 
