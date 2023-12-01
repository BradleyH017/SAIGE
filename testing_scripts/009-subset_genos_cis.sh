#!/usr/bin/env bash
#### Bradley 2023
#### Performing genome-wide analysis of trans-effects from cis-variants using SAIGE

# Load modules and docker
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
condition_col=""
condition=""
covariates="age_imputed,sex,Keras:predicted_celltype_probability"
covariates_cell="Keras:predicted_celltype_probability"
expression_pca=True
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5
trans_by_cis=true
cis_nominal_p_cut_off="5e-8"

# Set up dir
if [ -n "$condition_col" ]; then
    catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
    catdir=${general_file_dir}/${aggregate_on}/${level}
fi 

# First need to subset the genotype files for those that pass significance (may also want to do tests for independence before this)
for file in "${catdir}"/*.txt.gz; do
    # Extract the file name without the extension
    filename=$(basename "$file" .txt.gz)
    # Use zcat to read the compressed file, awk to filter rows, and grep to select columns
    zcat "$file" | awk -v col="13" -v val="$cis_nominal_p_cut_off" '$col < val' | cut -f "3" > "${filename}_cis_below_threshold.txt"
done

# Then want to combine
for chr_num in {1..22}; do
    echo $chr_num
    cat ${catdir}/chr${chr_num}_nPC_*_cis_below_threshold.txt >> ${catdir}/all_cis_below_threshold.txt
done

# Then subset the genotypes for these variants (bcf conda)
geno_dir=${general_file_dir}/genotypes
plink --bfile ${geno_dir}/plink_genotypes --extract ${catdir}/all_cis_below_threshold.txt --make-bed --out ${geno_dir}/plink_genotypes_cis_${level}


