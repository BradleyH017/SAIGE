#!/usr/bin/env bash
#### Bradley 2023
#### Summary of the iteration of results across chromosome 1

# Define options for this test (will ultimately be inherited) and general options
level="Enterocyte"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="category__machine"
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

# Set up dir
if [ -n "$condition_col" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

# Concat all of the ACAT results from chromosome 1 within each nPC
for n_expr_pcs in 5 10 15 20 $knee; do
        echo $n_expr_pcs
        for f in ${catdir}/*_npc${n_expr_pcs}_cis_ACAT.txt
        do
                tail -n +2 $f >> ${catdir}/chr1_nPC_${n_expr_pcs}_ACAT_all.txt
        done
done

# Execute R script to make summary plots [TO DO] and decide which nPC value to use
# scvi_gpu3 env
Rscript testing_scripts/bin/sum_ngeno_pcs_iteration_chrom1.R -c $catdir

# Load optimum PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
optim_npcs=$(<"$optim_npcs_file")

# Remove files for which the non-optimimum nPC is used
# Iterate over the desired values
for n_expr_pcs in 5 10 15 20 $knee; do
        echo $n_expr_pcs
        if [ "$n_expr_pcs" -ne "$optim_npcs" ]; then
                # Remove files that don't match the current value of n_optim_pcs
                rm ${catdir}/*npc${n_expr_pcs}*
                rm ${catdir}/*nPC_${n_expr_pcs}*
        fi
done

# Combine the raw results for this nPC
for f in ${catdir}/*_npc${optim_npcs}_cis.txt
        do
        tail -n +2 $f >> ${catdir}/chr1_nPC_${optim_npcs}.txt
done

# Remove intermediate files
rm ${catdir}/ENSG*

