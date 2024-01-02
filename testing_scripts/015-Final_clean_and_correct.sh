#!/usr/bin/env bash
#### Bradley 2023
#### Cleaning and correcting p-values from the trans-by-cis analysis (across chromosome)
# bsub -o logs/saige_trans_final-%J-output.log -e logs/saige_trans_final-%J-error.log -q long -G team152 -n 1 -M 40000 -a "memlimit=True" -R "select[mem>40000] rusage[mem=40000] span[hosts=1]" -J "saige_trans_final" < testing_scripts/015-Final_clean_and_correct.sh 

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

echo "~~~~~~~~~~~~~~~~~~~~~~~~ Working on:${level} ~~~~~~~~~~~~~~~~~~~~~~~~"

# Load n_expr_pcs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Delete the unwanted trans-eQTL files
# rm ${catdir}/trans/ENSG*
for c in {1..22}; do
    echo $c
    rm ${catdir}/trans/chr${c}_nPC_${n_expr_pcs}_trans_by_cis.txt # This one has not been filtered for the trans-only effects
done

# Remove the per-gene step1 files and other unwanted files
echo "Removing temporary files"
find ${catdir} -type f -name '*_region*' -print0 | xargs -0 rm
find ${catdir} -type f -name '*varianceRatio*' -print0 | xargs -0 rm
find ${catdir} -type f -name '*.rda' -print0 | xargs -0 rm
find ${catdir} -type f -name '*_cis_ACAT.txt' -print0 | xargs -0 rm
find ${catdir} -type f -name '*_cis_minimum_q.txt' -print0 | xargs -0 rm
find ${catdir} -type f -name 'ENSG*' -print0 | xargs -0 rm
find ${catdir} -type f -name '*all_genes_significant_genes.txt' -print0 | xargs -0 rm

# Remove the per chromosome list of genes
for c in {1..22} X Y MT; do
    rm ${catdir}/chr${c}_genes.txt
done
rm ${catdir}/non_chr1_genes.txt

# Move the other input files into a miscellaneous directory
mkdir -p ${catdir}/run_params
for f in ${catdir}/knee.txt ${catdir}/optim_nPCs_chr1.txt ${catdir}/test_genes.txt ${catdir}/expr_nPCs_check_chr1.txt; do
    mv $f ${catdir}/run_params
done

# Reheader the other cis files
echo "Reheading the cis results:"
for c in {1..22}; do
    echo "Reheading results for chromsome ${c}"
    # ACAT - all
    echo -e "Gene\tACAT_p\ttop_Marker_ID\ttop_p" >> ${catdir}/${c}_test.tsv
    cat ${catdir}/chr${c}_nPC_${n_expr_pcs}_ACAT_all.txt >> ${catdir}/${c}_test.tsv && mv ${catdir}/${c}_test.tsv ${catdir}/chr${c}_nPC_${n_expr_pcs}_ACAT_all.txt
    # minimum q all
    echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tp.value.NA\tIs.SPA\tN\tGene\tqvalues_within_genes\tlfdr_within_genes\tqvalues_across_genes\tlfdr_across_genes" >> ${catdir}/${c}_test.tsv
    tail -n +2 ${catdir}/chr${c}_nPC_${n_expr_pcs}_minimum_q_all_genes.txt >> ${catdir}/${c}_test.tsv && mv ${catdir}/${c}_test.tsv  ${catdir}/chr${c}_nPC_${n_expr_pcs}_minimum_q_all_genes.txt
    # all cis results
    echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tp.value.NA\tIs.SPA\tN\tGene\tqvalues_within_genes\tlfdr_within_genes" >> ${catdir}/${c}_test.tsv
    cat ${catdir}/chr${c}_nPC_${n_expr_pcs}.txt >> ${catdir}/${c}_test.tsv && mv ${catdir}/${c}_test.tsv ${catdir}/chr${c}_nPC_${n_expr_pcs}.txt
done

# Make a directory for these and move them
echo "Moving cis files"
mkdir -p ${catdir}/cis
mv ${catdir}/*.txt ${catdir}/cis

# Also compress the per-chromosome summary stats
for c in {1..22}; do
    echo "Compressing results for chromsome ${c}"
    gzip ${catdir}/cis/chr${c}_nPC_${n_expr_pcs}.txt
done

# Re-header the conditional file
echo "Reheading the conditional file"
echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tp.value.NA\tIs.SPA\tBETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\tp.value.NA_c\tN\tqvalues_c\tlfdr_c\tround\tGene" >> ${catdir}/conditional/temp.tsv
cat ${catdir}/conditional/all_conditionally_independent_effects.txt >> ${catdir}/conditional/temp.tsv && mv ${catdir}/conditional/temp.tsv ${catdir}/conditional/all_conditionally_independent_effects.txt

# Remove unwanted conditional files
echo "Removing temporary files"
rm ${catdir}/conditional/ENSG*
rm ${catdir}/conditional/chr*

# Remove the input files
rm -r ${catdir}/per_gene_input_files

# Compress the trans-eQTL summary stats
for c in {1..22}; do
    echo "Compressing trans results for chromsome ${c}"
    gzip ${catdir}/trans/chr${c}_nPC_${n_expr_pcs}_trans_by_cis_no_cis.txt
done

# Tidy up the genotypes
echo "Removing the unneccessary genotype file"
rm ${general_file_dir}/genotypes/plink_genotypes_cis_${level}*

# Tidy up the last remaining files
echo "Removing last remaining files"
rm ${catdir}/*size_temp
rm ${catdir}/trans/*__npc*
rm ${catdir}/cis/*_temp_independent.txt