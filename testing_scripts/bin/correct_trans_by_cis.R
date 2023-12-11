#!/usr/bin/env Rscript

###### Bradley Dec 2023
library(optparse)
library(dplyr)

# Inherit options
option_list = list(
    make_option(c("-c", "--catdir"), action = "store", default = NA, type ="character",
        help="Where are the files located?"),
    make_option(c("-chr", "--chromosome"), action = "store", default = NA, type ="character",
        help="Which chromosome to work on?"),
    make_option(c("-n", "--n_expr_pcs"), action = "store", default = NA, type ="character",
        help="What is the optimum number of expression PCs to use"),
    make_option(c("-w", "--cis_window"), action = "store", default = NA, type ="character",
        help="What is the window to define cis-eQTL vs trans-eQTL"),
    make_option(c("-a", "--annotation_file"), action = "store", default = NA, type ="character",
        help="Annotation file used for the analysis")
    )
opt = parse_args(OptionParser(option_list = option_list))
catdir = opt$c
chr = opt$chr
n_expr_pcs = opt$n
window = opt$w
annotation_file = opt$a

# Testing: 
# catdir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/label__machine/Tuft_cell"
# chr = "1"

# Load in the results for this chromosome
trans = read.delim(sprinf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis.txt"), chr, n_expr_pcs), header=F)
colnames(trans)=c("CHR", "POS",  "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA", "Is.SPA", "N", "Gene")

# Check the position of these variants compared with the trans
gtf=read.delim(annotation_file)
# Merge these
gtf = gtf[gtf$feature_id %in% trans$Gene,]
colnames(gtf)[2:ncol(gtf)] = paste0("gene_", colnames(gtf)[2:ncol(gtf)]) 
colnames(gtf)[1] = "Gene"
trans = merge(trans, gtf, by="Gene")
# Subset cis-effects
filter_rows <- function(row) {
  chromosome_match <- row["CHR"] == row["gene_chromosome"]
  position_match <- abs(row["POS"] - row["gene_start"]) < window
  return(!(chromosome_match & position_match))
}
filtered_trans <- trans[apply(trans, 1, filter_rows), ]
# Save the filtered output
write.table(filtered_trans, sprinf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis.txt"), chr, n_expr_pcs))

# Perform BH correction within genes (ISSUE: THIS ASSUMES INPUT VARIANTS ARE INDEPENDENT - WHICH AS OF YET IS NOT TRUE)
filtered_trans <- filtered_trans %>%
  group_by(Gene) %>%
  mutate(p.value.bonf = p.adjust(p.value, method = "bonferroni"))

trans_bonf = filtered_trans %>%
  group_by(Gene) %>%
  slice_min(order_by = p.value.bonf)

trans_bonf$FDR = p.adjust(trans_bonf$p.value.bonf, method="fdr")

# Save the output
write.table(trans_bonf, sprinf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis_bonf_fdr.txt"), chr, n_expr_pcs))
