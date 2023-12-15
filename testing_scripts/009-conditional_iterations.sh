#!/usr/bin/env bash
# Perform the SAIGEQTL analysis of single cell expression from TI - conditioning on the top variant (present in the ACAT test results) - run per test gene
# bsub -o logs/saige_conditional-%J-%I-output.log -e logs/saige_conditional-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_conditional[1-12899]%500" < testing_scripts/009-conditional_iterations.sh 

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

# If directory does not exist yet, make it
cond_dir=${catdir}/conditional
mkdir -p $cond_dir

# Find the minimum q-value for this gene
top_q=$(awk -F'\t' 'NR==2 {print $18}' ${step2prefix}_minimum_q.txt)

# If this is above 0.05, finish script here
threshold=0.05
if (( $(echo "$top_q > $threshold" | bc -l) )); then
    echo "Not performing conditional analysis: q-value for first pass > $threshold"
    exit 0  # 0 means success in Bash
fi

# If not, perform conditional
topvariant=$(awk -F'\t' 'NR==2 {print $3}' ${step2prefix}_minimum_q.txt)

# Re-derive the gene_chr
gene_chr=$(awk -F'\t' 'NR==2 {print $1}' ${step2prefix}_minimum_q.txt)

# Perform conditional cis-eQTL analysis on this variant, for this gene
# Output is called 'round2' as the first test was round 1
echo "~~~~~~~~~~~~~~~~PERFORMING THE SECOND ROUND OF EQTL ANALYSIS~~~~~~~~~~~~~~~~~"
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
        --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
        --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
        --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
        --SAIGEOutputFile=${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round2.txt    \
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

# Perform qvalue correction of these results
Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round2.txt -c "20" -n "c-${topvariant}.qvalues" -w "TRUE"

# Now perform an additional 3 rounds if we keep finding conditional 
for c in {3..5}; do
    echo "~~~~~~~~~~~~~~~~PERFORMING ROUND ${c} OF EQTL ANALYSIS~~~~~~~~~~~~~~~~~"
    previous=$((c - 1))
    # Find the minimum q of the round previously and whether this passes the threshold [new conditional column is always column 23]
    top_q_row=$(tail -n +2 "${catdir}/conditional/${gene}__npc${n_expr_pcs}_cis_round${previous}.txt" | sort -t$'\t' -k23,23n | head -n 1)
    top_q=$(echo "$top_q_row" | awk -F' ' '{print $23}')
    # Stop if no further independent effects 
    if (( $(echo "$top_q > $threshold" | bc -l) )); then
        echo "Not performing conditional analysis round ${c}: q-value for last pass > $threshold"
        exit 0  # 0 means success in Bash
    fi 
    # But if not, combine the new top variant with the previous set (3rd column)
    new_top=$(echo "$top_q_row" | awk -F' ' '{print $3}')
    # Combine to condition on both
    topvariant="${topvariant},${new_top}"
    # Perform SAIGE, now conditioning on these two variants
    singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
        --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
        --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
        --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
        --SAIGEOutputFile=${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round${c}.txt    \
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
    # Perform correction of results
    Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${cond_dir}/${gene}__npc${n_expr_pcs}_cis_round${c}.txt -c "20" -n "c-${topvariant}.qvalues" -w "TRUE"
    echo "~~~~~~~~~~~~~~~~DONE CONDITIONAL ROUND ${c}!~~~~~~~~~~~~~~~~"
done
