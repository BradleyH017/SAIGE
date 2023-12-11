#!/usr/bin/env bash
# Perform the SAIGEQTL analysis of single cell expression from TI - conditioning on the top variant (present in the ACAT test results) - run per test gene
# bsub -o logs/saige_array_test-%J-%I-output.log -e logs/saige_array_test-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_array_test[1-1542]%300" < testing_scripts/005-run_SAIGE_1_2_3_chrom1.sh 

# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
level="Goblet_cell_middle_villus"
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

# Define the number of expression PCs to use (from the optimum)
if [[ "$expression_pca" == "true" ]]; then
    # Load the optimum number of PCs (expression)
    optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
    n_expr_pcs=$(<"$optim_npcs_file")
else
    n_expr_pcs=0
fi

# Load the gene to test
gene_list=${catdir}/test_genes.txt
gene=$(head $gene_list -n ${LSB_JOBINDEX} | tail -n 1)

# Define the prefixes
step1prefix=${catdir}/${gene}_npc${n_expr_pcs}
step2prefix=${catdir}/${gene}__npc${n_expr_pcs}_cis

for round in {2..5}; do

# Find the top variant per test gene
topvariant=

# Perform conditional cis-eQTL analysis on this variant, for this gene
# Output is called 'round2' as the first test was round 1
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
        --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
        --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
        --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
        --SAIGEOutputFile=${step2prefix}_round2.txt     \
        --chrom=${gene_chr}       \
        --minMAF=0.05 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --is_imputed_data=TRUE \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${step2prefix}_region_file.txt     \
        --markers_per_chunk=10000 \
        --condition=$topvariant

# The order of variants from this are the same as those in the original cis-results. So can just paste the interesting columns (16:22) onto the original results file for this gene
suffix="_round2_${topvariant}"
awk -F'\t' -v OFS='\t' -v suffix="$suffix" 'NR == 1 {for (i=16; i<=22; i++) $i = $i suffix} 1' "${step2prefix}_round2.txt" > ${gene}_round2_to_add.txt
paste <(awk -F'\t' 'NR==1 {for (i=1; i<=22; i++) print $i}' "${step2prefix}.txt") <(awk -F'\t' 'NR > 1 {print $16, $17, $18, $19, $20, $21, $22}' ${gene}_round2_to_add.txt) > test.txt


# Select the variant with the minimum conditional p-value, append this to the $topvariant variable and repeat






                # Add gene name to the output
                awk -v new_val=${gene} 'BEGIN {OFS="\t"} {print $0, new_val}' "${step2prefix}.txt" > tmp && mv tmp ${step2prefix}.txt

awk -F'\t' -v col="3" '{print $col}' "${step2prefix}_round2.txt"
