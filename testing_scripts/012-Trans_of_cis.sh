#!/usr/bin/env bash
#### Bradley 2023
#### Performing genome-wide analysis of trans-effects from cis-variants using SAIGE
# bsub -o logs/saige_trans-%J-%I-output.log -e logs/saige_trans-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_trans[1-12899]%500" < testing_scripts/012-Trans_of_cis.sh 


# Load modules and docker
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


echo "Prepping the directory variables"
# Construct the level directory path
if [ -n "$condition_col" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

# If submitting as an array job for all genes:
gene_list=${catdir}/test_genes.txt
gene=$(head $gene_list -n ${LSB_JOBINDEX} | tail -n 1)

# Add genotype PCs to the covariates
for ((i=1; i<=n_geno_pcs; i++)); do
        covariates+=",PC$i"
done

# Define the number of expression PCs to use
if [[ "$expression_pca" == "true" ]]; then
    # Load the optimum number of PCs (expression)
    optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
    n_expr_pcs=$(<"$optim_npcs_file")
else
    n_expr_pcs=0
fi

# Execute SAIGE for this value of expression PCs
echo "Working for expression PC number"
echo $n_expr_pcs
echo "Prepping the covariates"

 # Generate the expression PC' string up to the numeric value
if [ "${n_expr_pcs}" != 0 ]; then
    covariates_cell_pc="${covariates_cell},$(printf "xPC%d," $(seq "$n_expr_pcs") | sed 's/,$//')"
else
    covariates_cell_pc="$covariates_cell"
fi

# Fix covariate issue (replacing ':' with '_') in both the covariates and covatiates_cell
covariates="${covariates//:/_}"
covariates_cell_pc="${covariates_cell_pc//:/_}"
# Combine
covariates_sample_cell=$(echo "$covariates,$covariates_cell_pc" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
# Specify sample-level covariates
covariates_sample=$(echo "$covariates" | tr ',' '\n' | sort | comm -23 - <(echo "$covariates_cell_pc" | tr ',' '\n' | sort) | tr '\n' ',' | sed 's/,$//')

# Make trans dir if not present
trans_dir=${catdir}/trans
mkdir -p $trans_dir

# NOw perform the cis-eqtl-wide test
echo "Testing eQTLs"
step1prefix=${catdir}/${gene}_npc${n_expr_pcs}
step2prefix=${trans_dir}/${gene}__npc${n_expr_pcs}_trans_by_cis

# Now perform the Genome-wide analysis - NOTE: This also includes cis-variants
echo "Performing the trans-by-cis-only analysis"
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
    --bedFile=${general_file_dir}/genotypes/plink_genotypes_cis_${level}.bed      \
    --bimFile=${general_file_dir}/genotypes/plink_genotypes_cis_${level}.bim      \
    --famFile=${general_file_dir}/genotypes/plink_genotypes_cis_${level}.fam      \
    --SAIGEOutputFile=${step2prefix}.txt     \
    --minMAF=0.05 \
    --minMAC=20 \
    --LOCO=FALSE    \
    --GMMATmodelFile=${step1prefix}.rda     \
    --SPAcutoff=2 \
    --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
    --markers_per_chunk=10000

# Add gene name to the output
awk -v new_val=${gene} 'BEGIN {OFS="\t"} {print $0, new_val}' "${step2prefix}.txt" > tmp_${gene} && mv tmp_${gene} ${step2prefix}.txt 

echo "Finished analysis, removed intermediate files"