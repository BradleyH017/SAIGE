#!/usr/bin/env bash
#### Bradley 2023
#### Cleaning and correcting p-values from the trans-by-cis analysis (per chromosome)
# bsub -o logs/saige_tidy_trans-%J-%I-output.log -e logs/saige_tidy_trans-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_tidy_trans[1-22]" < testing_scripts/011-Clean_trans_of_cis_correct.sh 


module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
level="Tuft_cell"
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
expression_pca=True
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5

# Construct the category directory path
if [ "$condition_col" != "NULL" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

#  Use the job ID as the chromosome number in this instance
gene_chr=${LSB_JOBINDEX}

# Group all results from the same chromosome into the same results file
# Load the optimum number of PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Tidy up the results from this chromosome
echo "Tidying up the output"
while read gene; do
        # Gather the trans-by-cis results
        tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_trans_by_cis.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_trans_by_cis.txt
        # Gather the trans-by-cis ACAT results
        tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_trans_by_cis_ACAT.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_trans_by_cis_ACAT_all.txt
        # remove the intermediate files
        rm ${catdir}/${gene}*
done < ${catdir}/chr${gene_chr}_genes.txt

# Perform correction of the trans-by-cis-results
echo "Correcting the trans-results"
Rscript testing_scripts/bin/correct_trans_by_cis.R -c $catdir -chr $gene_chr -n $n_expr_pcs -w $cis_window -a $annotation__file