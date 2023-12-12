#!/usr/bin/env bash
level="T_cell_CD8_1"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="label__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=10
condition_col="NULL" #Specify 'NULL' if want to include all cells
condition="NULL"
covariates="age_imputed,sex,Keras:predicted_celltype_probability"
covariates_cell="Keras:predicted_celltype_probability"
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
scale_covariates=false # SAIGE has an option to do this - leave it to SAIGE
cis_window=1000000

# Run script to subset anndata object based on this aggregation column, identify genes to test, make the neccessary input files for SAIGE
# bsub -o logs/saige_prep-%J-output.log -e logs/saige_prep-%J-error.log -q long -G team152 -n 4 -M 70000 -a "memlimit=True" -R "select[mem>70000] rusage[mem=70000] span[hosts=1]" -J "saige_prep" < testing_scripts/003-prep_adata_per_resolution.sh 
# mamba activate scvi_gpu3
python testing_scripts/bin/prep_adata_saige.py \
    --phenotype__file $phenotype__file \
    --aggregate_on $aggregate_on \
    --genotype_pc__file $genotype_pc__file \
    --genotype_id $genotype_id \
    --sample_id $sample_id \
    --general_file_dir $general_file_dir \
    --nperc $nperc \
    --condition_col $condition_col \
    --condition $condition_col \
    --covariates $covariates \
    --scale_covariates $scale_covariates \
    --expression_pca $expression_pca \
    --level $level