#!/usr/bin/env bash
#### Bradley 2023
#### Cleaning and correcting p-values from the trans-by-cis analysis (across chromosome)
# bsub -o logs/saige_trans_final-%J-output.log -e logs/saige_trans_final-%J-error.log -q long -G team152 -n 1 -M 40000 -a "memlimit=True" -R "select[mem>40000] rusage[mem=40000] span[hosts=1]" -J "saige_trans_final" < testing_scripts/014-trans_gather_correct_clean.sh 


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

# Load the optimum number of PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Gather these results
echo "Gathering results for:"
head -n 1 ${catdir}/trans/chr1_nPC_${n_expr_pcs}_trans_by_cis_no_cis_bonf.txt >> ${catdir}/trans/all_nPC_${n_expr_pcs}_trans_by_cis_no_cis_bonf.txt
for c in {1..22}; do
    echo $c
    tail -n +2 ${catdir}/trans/chr${c}_nPC_${n_expr_pcs}_trans_by_cis_no_cis_bonf.txt >> ${catdir}/trans/all_nPC_${n_expr_pcs}_trans_by_cis_no_cis_bonf.txt
done

# Then perform qvalue correction of this file and overwrite
echo "Performing qvalue correction across all genes"
Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${catdir}/trans/all_nPC_${n_expr_pcs}_trans_by_cis_no_cis_bonf.txt -c "22" -n "bonf_qvalue_across_genes" -w "FALSE"

# DONE!