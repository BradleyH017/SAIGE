#!/usr/bin/env bash
# Cleaning up the conditional analysis results and combining into per-chromosome results
# bsub -o logs/saige_tidy_conditional-%J-%I-output.log -e logs/saige_tidy_conditional-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_tidy_conditional[1-22]" < testing_scripts/010-clean_combine_conditional.sh 


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

# Run per chromosome within category
gene_chr=${LSB_JOBINDEX}

# Across each gene from this chromosome with significant effects
while read gene; do
    echo $gene
    # Across each gene of this chromosome
    # Grab the round 1 results
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"","\t\t\t\t",$16,$18,$19}' ${catdir}/${gene}__npc${n_expr_pcs}_cis_minimum_q.txt | tail -n 1 >> ${catdir}/${gene}_temp_independent.txt
    # For every round for this gene in the conditional analysis directory (this will automatically be in order of round)
    for f in ${catdir}/conditional/${gene}*_minimum_q*; do
        tail -n 1 $f >> ${catdir}/${gene}_temp_independent.txt
    done
    # Add the 'round' column to the output
    awk 'BEGIN {OFS="\t"} {print $0, NR}' ${catdir}/${gene}_temp_independent.txt > ${catdir}/${gene}_temp_independent2.txt && mv ${catdir}/${gene}_temp_independent2.txt ${catdir}/${gene}_temp_independent.txt
    # Add the gene name to the file before combining with other genes
    awk -v var_value="$gene" 'BEGIN {OFS="\t"} {print $0, var_value}' ${catdir}/${gene}_temp_independent.txt > ${catdir}/${gene}_temp_independent2.txt && mv ${catdir}/${gene}_temp_independent2.txt ${catdir}/${gene}_temp_independent.txt
    # Append this onto those for the whole chromsome 
    head -n 10 ${catdir}/${gene}_temp_independent.txt >> ${catdir}/conditional/chr${gene_chr}_conditionally_independent_effects.txt 
    # Remove all of the conditional files for this gene and the temporary files
    #rm ${catdir}/conditional/${gene}*
    #rm ${catdir}/${gene}_temp_independent.txt
done <${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_minimum_q_all_genes_significant_genes.txt