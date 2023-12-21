#!/usr/bin/env bash
#### Bradley 2023
#### Performing genome-wide analysis of trans-effects from cis-variants using SAIGE

# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
level="Pericytes"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="label__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=1
condition_col=""
condition=""
covariates="age_imputed,sex"
covariates_cell=""
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5
use_GRM=FALSE
repo_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE"


# Set up dir
if [ -n "$condition_col" ]; then
    catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
    catdir=${general_file_dir}/${aggregate_on}/${level}
fi 

# Load optimum PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Combine the conditional results across chromosomes into a single conditional all file
for c in {1..22}; do
    echo $c
    cat ${catdir}/conditional/chr${c}_conditionally_independent_effects.txt >> ${catdir}/conditional/all_conditionally_independent_effects.txt
done

# Subset for those that are qvalue < 0.05 (column 23)
awk -v col="23" -v threshold="0.05" '$col < threshold' ${catdir}/conditional/all_conditionally_independent_effects.txt > ${catdir}/conditional/all_conditionally_independent_effects_q_less_0.05.txt

# Extract just variants and save to file
awk '{print $3}' ${catdir}/conditional/all_conditionally_independent_effects_q_less_0.05.txt > ${catdir}/conditional/all_conditionally_independent_effects_q_less_0.05_variants.txt

# Then subset the genotypes for these variants (bcf conda)
geno_dir=${general_file_dir}/genotypes
plink --bfile ${geno_dir}/plink_genotypes_ordered --extract ${catdir}/conditional/all_conditionally_independent_effects_q_less_0.05_variants.txt --make-bed --out ${geno_dir}/plink_genotypes_cis_${level}


