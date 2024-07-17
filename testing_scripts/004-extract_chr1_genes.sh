#!/usr/bin/env bash
# Extracting genes on chromosome one to allow iteration over these first
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
covariates="age_imputed,sex,total_counts"
covariates_cell="total_counts"
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000

echo "Prepping the directory variables"
# Construct the category directory path
if [ "$condition_col" != "NULL" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

# Extract chromosome 1 genes
while read gene; do
    gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
    if [ "$gene_chr" == "1" ]; then
        echo "$gene" >> ${catdir}/chr1_genes.txt
    else
        echo "$gene" >> ${catdir}/non_chr1_genes.txt
    fi
done <${catdir}/test_genes.txt

# For clearing up later (also make lists for non-chromosome 1 genes per chromosome)
while read gene; do
    gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
    echo "$gene" >> ${catdir}/chr${gene_chr}_genes.txt
done <${catdir}/non_chr1_genes.txt