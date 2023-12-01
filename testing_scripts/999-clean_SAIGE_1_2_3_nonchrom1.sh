#!/usr/bin/env bash
#### Bradley 2023
# Cleaning the ouput of SAIGE on the rest of the chromosomes

# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
category="Tuft_cell"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="label__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=20
condition_col=""
condition=""
covariates="age_imputed,sex,Keras:predicted_celltype_probability"
covariates_cell="Keras:predicted_celltype_probability"
expression_pca=True
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000

# Set up dir
if [ -n "$condition_col" ]; then
    catdir=${general_file_dir}/${aggregate_on}/${category}/${condition_col}/${condition}
else
    catdir=${general_file_dir}/${aggregate_on}/${category}
fi

# Load optimum PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_geno_pcs=$(<"$optim_npcs_file")

# Gather results for each chromosome (also compressing these)
echo "Cleaning the output of the rest of the chromosomes"
while read gene; do
    echo $gene
    gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
    # Gather the cis results
    tail -n +2 ${catdir}/${gene}__npc${n_geno_pcs}_cis.txt >> ${catdir}/chr${gene_chr}_nPC_${n_geno_pcs}.txt
    # Gather the ACAT results
    tail -n +2 ${catdir}/${gene}__npc${n_geno_pcs}_cis_ACAT.txt >> ${catdir}/chr${gene_chr}_nPC_${n_geno_pcs}_ACAT_all.txt
    # remove the intermediate files
    rm ${catdir}/${gene}*
done <${catdir}/non_chr1_genes.txt

# Gzip each chromosome
for chr_num in {1..22}; do
    gzip ${catdir}/chr${chr_num}_nPC_${n_geno_pcs}.txt
done


