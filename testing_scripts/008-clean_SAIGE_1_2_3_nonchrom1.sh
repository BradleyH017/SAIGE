#!/usr/bin/env bash
#### Bradley 2023
#### Cleaning SAIGE on the rest of the chromomes (per chromosome) - This needs to be submitted as an array - 1 job per non '1' chromosome
# bsub -o logs/saige_tidy-%J-%I-output.log -e logs/saige_tidy-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_tidy[1-21]" < testing_scripts/008-clean_SAIGE_1_2_3_nonchrom1.sh 


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
covariates="age_imputed,sex,Keras:predicted_celltype_probability"
covariates_cell="Keras:predicted_celltype_probability"
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5
repo_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE"

# Construct the category directory path
if [ "$condition_col" != "NULL" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

#  Use the job ID as the chromosome number in this instance (running on non chromosome 1 so add 1)
gene_chr=$((LSB_JOBINDEX + 1))

# Load the optimum number of PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Tidy up the results from this chromosome
while read gene; do
        # Gather the minimum q results
        tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_cis_minimum_q.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_minimum_q_all_genes.txt
        # Gather the cis results
        tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_cis.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}.txt
        # Gather the ACAT results
        tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_cis_ACAT.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_ACAT_all.txt
        # remove the intermediate files
        #rm ${catdir}/${gene}*
done < ${catdir}/chr${gene_chr}_genes.txt

# Then calculate the q-value for these results 
Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_minimum_q_all_genes.txt -c "18" -n "qvalues_across_genes" -w "FALSE"

# Output list of Genes with q-value (across chromosome) significant effects
awk -v threshold="0.05" -F '\t' 'NR==1 || $20 < threshold {print $17}' ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_minimum_q_all_genes.txt | tr ',' '\t' | tail -n +2 > ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_minimum_q_all_genes_significant_genes.txt

# bsub -o logs/saige_tidy-%J-%I-output.log -e logs/saige_tidy-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_tidy[1-21]" < testing_scripts/008-clean_SAIGE_1_2_3_nonchrom1.sh 
