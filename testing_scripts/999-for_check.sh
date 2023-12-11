echo "Estimating the variance"
use_GRM=FALSE
# saige_filt_expr_input.txt = expression and covariate matrix per cell
# catdir = results directory
# n_expr_pcs = number of expression PCs we are running here
# gene = gene to run on - present within the header of saige_filt_expr_input.txt
# annotation_file = gtf file to grab gene coordinates from
# cis_window = distance from TSS to run cis-analysis within
# covariates_sample = The sample covariates to include in the analysis = age_imputed,PC1,PC2,PC3,PC4,PC5,sex
# covariates_sample_cell = Combination of both sample- and cell-level covariates included in analysis [xPC = expression PC, PC = genotype PC]= age_imputed,PC1,PC2,PC3,PC4,PC5,sex,xPC1,xPC10,xPC2,xPC3,xPC4,xPC5,xPC6,xPC7,xPC8,xPC9
singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R \
    --useSparseGRMtoFitNULL=FALSE  \
    --useGRMtoFitNULL=$use_GRM \
    --phenoFile=${catdir}/saige_filt_expr_input.txt	\
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
    # Find the coordinates/chromosome of the given gene
    gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
    gene_start=$(awk -v search="$gene" '$1 == search {print $2}' "$annotation__file")
    end_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {print value + add}')
    start_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {new_value = value - add; if (new_value < 0) new_value = 0; print new_value}')
    # Save this into a temporary file to perform the cis-only analysis
    echo -e "${gene_chr}\t${start_region}\t${end_region}" > ${step2prefix}_region_file.txt

    # Now perform the cis-only analysisq2
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