#!/usr/bin/env bash
# Perform the SAIGEQTL analysis of single cell expression from TI (test)
# bsub -o logs/saige_array_test-%J-%I-output.log -e logs/saige_array_test-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "saige_array_test[1-1305]%400" < testing_scripts/005-run_SAIGE_1_2_3_chrom1.sh 

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

# If submitting as an array job for chromosome 1 genes only:
gene_list=${catdir}/chr1_genes.txt
gene=$(head $gene_list -n ${LSB_JOBINDEX} | tail -n 1)

# Add genotype PCs to the covariates
for ((i=1; i<=n_geno_pcs; i++)); do
        covariates+=",PC$i"
done

# Perform SAIGE with default options on a number of genes from the gene list file - looping over the number of expression PCs from 5,10,15, 20 and the knee if using expression PCs
# Load knee
knee_file=${catdir}/knee.txt
knee=$(<"$knee_file")

if [[ "$expression_pca" == "true" ]]; then
        n_expr_pc_params="0 5 10 15 20 $knee"
else
        n_expr_pc_params=0
fi

for n_expr_pcs in $n_expr_pc_params; do
        echo "Working for expr PC number"
        echo $n_expr_pcs
        echo "Prepping the covariates"

        # Generate the expression PC' string up to the numeric value
        if [ "${n_expr_pcs}" != 0 ]; then
                covariates_cell_pc="${covariates_cell},$(printf "xPC%d," $(seq "$n_expr_pcs") | sed 's/,$//')"
        else
                covariates_cell_pc="$covariates_cell"
        fi

        # Fix covariate issue (replacing ':' with '_') in both the covariates and covatiates_cell
        covariates="${covariates//:/_}"
        covariates_cell_pc="${covariates_cell_pc//:/_}"
        # Combine
        covariates_sample_cell=$(echo "$covariates,$covariates_cell_pc" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
        # Specify sample-level covariates
        covariates_sample=$(echo "$covariates" | tr ',' '\n' | sort | comm -23 - <(echo "$covariates_cell_pc" | tr ',' '\n' | sort) | tr '\n' ',' | sed 's/,$//')

        # Fix input of first ','
        if [[ $covariates_sample_cell == ,* ]]; then
        # Remove the first character
                covariates_sample_cell=${covariates_sample_cell:1}
        fi

        echo "Estimating the variance"
        singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R \
                --useSparseGRMtoFitNULL=FALSE  \
                --useGRMtoFitNULL=$use_GRM \
                --phenoFile=${catdir}/per_gene_input_files/${gene}_saige_filt_expr_input.txt	\
                --phenoCol=${gene}       \
                --covarColList=${covariates_sample_cell}    \
                --sampleCovarColList=${covariates_sample}      \
                --sampleIDColinphenoFile=${genotype_id} \
                --traitType=count \
                --outputPrefix=${catdir}/${gene}_npc${n_expr_pcs} \
                --skipVarianceRatioEstimation=FALSE  \
                --isRemoveZerosinPheno=FALSE \
                --isCovariateOffset=FALSE  \
                --isCovariateTransform=TRUE  \
                --skipModelFitting=FALSE  \
                --tol=0.00001   \
                --sexCol="sex" \
                --FemaleCode="0" \
                --MaleCode="1" \
                --plinkFile=${general_file_dir}/genotypes/plink_genotypes_ordered      \
                --IsOverwriteVarianceRatioFile=TRUE

        # Perform the analysis cis-only or genome-wide
        echo "Testing eQTLs"
        step1prefix=${catdir}/${gene}_npc${n_expr_pcs}
        if [ "$cis_only" = true ]; then
                step2prefix=${catdir}/${gene}__npc${n_expr_pcs}_cis
                # Find the TSS/coordinates/chromosome of the given gene
                gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
                strand=$(awk -v search="$gene" '$1 == search {print $4}' "$annotation__file")
                # TSS is first number if strand is '+'
                if [ "$value" == "+" ]; then
                        gene_start=$(awk -v search="$gene" '$1 == search {print $2}' "$annotation__file")
                else
                        gene_start=$(awk -v search="$gene" '$1 == search {print $3}' "$annotation__file")
                fi
                end_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {print value + add}')
                start_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {new_value = value - add; if (new_value < 0) new_value = 0; print new_value}')
                # Save this into a temporary file to perform the cis-only analysis
                echo -e "${gene_chr}\t${start_region}\t${end_region}" > ${step2prefix}_region_file.txt

                # Now perform the cis-only analysis
                echo "Performing the cis-only analysis"
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
                        --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
                        --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
                        --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
                        --SAIGEOutputFile=${step2prefix}.txt     \
                        --chrom=${gene_chr}       \
                        --minMAF=0.05 \
                        --minMAC=20 \
                        --LOCO=FALSE    \
                        --is_imputed_data=TRUE \
                        --GMMATmodelFile=${step1prefix}.rda     \
                        --SPAcutoff=2 \
                        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
                        --rangestoIncludeFile=${step2prefix}_region_file.txt     \
                        --markers_per_chunk=10000
                
                # Also perform the ACAT test for this gene
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
                        --assocFile=${step2prefix}.txt      \
                        --geneName=$gene       \
                        --genePval_outputFile=${step2prefix}_ACAT.txt

                # Add gene name to the output
                awk -v new_val=${gene} 'BEGIN {OFS="\t"} {print $0, new_val}' "${step2prefix}.txt" > tmp && mv tmp ${step2prefix}.txt

                # Q-value correction of the per gene results (specify the input file, the column that encodes the p-value, the new column name and whether to run within gene)
                echo "Performing q-value correction"
                Rscript ${repo_dir}/testing_scripts/bin/qvalue_correction.R -f ${step2prefix}.txt -c "13" -n "qvalues" -w "TRUE"

                echo "Finished analysis"
        else
                step2prefix=${catdir}/${gene}__npc${n_expr_pcs}_gw.txt
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R &> ${catdir}/${gene}_npc${n_expr_pcs}_log.txt \
                        --bedFile=${general_file_dir}/genotypes/plink_genotypes.bed      \
                        --bimFile=${general_file_dir}/genotypes/plink_genotypes.bim      \
                        --famFile=${general_file_dir}/genotypes/plink_genotypes.fam      \
                        --SAIGEOutputFile=${step2prefix}     \
                        --minMAF=0.05 \
                        --minMAC=20 \
                        --LOCO=FALSE    \
                        --GMMATmodelFile=${step1prefix}.rda     \
                        --SPAcutoff=2 \
                        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
                        --markers_per_chunk=10000

                # Run the ACAT test on the genome wide results
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
                        --assocFile=${step2prefix}        \
                        --geneName=$gene       \
                        --genePval_outputFile=${catdir}/${gene}_gw_ACAT.txt
        fi
done        