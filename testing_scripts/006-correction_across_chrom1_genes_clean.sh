#!/usr/bin/env bash
# Performing multiple testing correction (across genes) for chromosome 1 and iterative conditional analysis for significant effects

# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
level="T_cell_CD4_CD40LGplus_2"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="label__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=10
condition_col=""
condition=""
covariates="age_imputed,sex"
covariates_cell=""
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000
n_geno_pcs=5
repo_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE"


echo "Prepping the directory variables"
# Construct the level directory path
if [ -n "$condition_col" ]; then
        catdir=${general_file_dir}/${aggregate_on}/${level}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${aggregate_on}/${level}
fi

# Load knee
knee_file=${catdir}/knee.txt
knee=$(<"$knee_file")

if [[ "$expression_pca" == "true" ]]; then
        n_expr_pc_params="0 5 10 15 20 $knee"
else
        n_expr_pc_params=0
fi

# Gather all of the qvalue results together (per nPC value)
for n_expr_pcs in $n_expr_pc_params; do
    echo "Grouping qvalue for nPC ${n_expr_pcs}"
    for f in ${catdir}/*_npc${n_expr_pcs}*cis_minimum_q*; do
        tail -n +2 ${f} >> ${catdir}/chr1_nPC_${n_expr_pcs}_minimum_q_all_genes.txt
    done 
    # Then calculate the q-value for these results (within chromosome 1)
    Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${catdir}/chr1_nPC_${n_expr_pcs}_minimum_q_all_genes.txt -c "18" -n "qvalues_across_genes" -w "FALSE"
done

# Then find the optimum number number of expression PCs to calculate the number of eGenes
Rscript ${repo_dir}/testing_scripts/bin/sum_ngeno_pcs_iteration_chrom1.R -c $catdir

# Load optimum PCs
optim_npcs_file=${catdir}/optim_nPCs_chr1.txt
n_expr_pcs=$(<"$optim_npcs_file")

# Gather results for this chromosome
echo "Cleaning the output of the rest of the chromosomes"
while read gene; do
    gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
    # Gather the cis results
    tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_cis.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}.txt
    # Gather the ACAT results
    tail -n +2 ${catdir}/${gene}__npc${n_expr_pcs}_cis_ACAT.txt >> ${catdir}/chr${gene_chr}_nPC_${n_expr_pcs}_ACAT_all.txt
done <${catdir}/chr1_genes.txt

# Remove the files without the optimum nPC
for n_expr_pcs_test in 0 5 10 15 20 $knee; do
        echo $n_expr_pcs_test
        if [ "$n_expr_pcs_test" -ne "$n_expr_pcs" ]; then
                # Remove files that don't match the current value of n_optim_pcs
                rm ${catdir}/*npc${n_expr_pcs_test}*
                rm ${catdir}/*nPC_${n_expr_pcs_test}*
        fi
done
