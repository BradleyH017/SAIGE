#!/usr/bin/env Rscript

###### Bradley Dec 2023
library(optparse)
library(dplyr)

# Inherit options
option_list = list(
    make_option(c("-d", "--catdir"), action = "store", default = NA, type ="character",
        help="Where are the files located?"),
    make_option(c("-c", "--chromosome"), action = "store", default = NA, type ="character",
        help="Which chromosome to work on?"),
    make_option(c("-n", "--n_expr_pcs"), action = "store", default = NA, type ="character",
        help="What is the optimum number of expression PCs to use"),
    make_option(c("-w", "--cis_window"), action = "store", default = NA, type ="character",
        help="What is the window to define cis-eQTL vs trans-eQTL"),
    make_option(c("-a", "--annotation_file"), action = "store", default = NA, type ="character",
        help="Annotation file used for the analysis")
    )
opt = parse_args(OptionParser(option_list = option_list))
catdir = opt$d
chr = opt$c
n_expr_pcs = opt$n
window = opt$w
annotation_file = opt$a

# Load in the results for this chromosome
print("reading in trans")
trans = read.delim(sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis.txt"), as.character(chr), as.character(n_expr_pcs)), header=F)
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
  position_match <- abs(as.numeric(row["POS"]) - as.numeric(row["gene_start"])) < window
  return(!(chromosome_match & position_match))
}
print("Filtering trans to remove out cis")
filtered_trans <- trans[apply(trans, 1, filter_rows), ]
# Save the filtered output
write.table(filtered_trans, sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis.txt"), as.character(chr), as.character(n_expr_pcs)))

# Perform BH correction within genes
print("Performing BH correction within gene")
filtered_trans <- filtered_trans %>%
  group_by(Gene) %>%
  mutate(p.value.bonf = p.adjust(p.value, method = "bonferroni"))

print("Performing FDR correction across genes")
trans_bonf = filtered_trans %>%
  group_by(Gene) %>%
  slice_min(order_by = p.value.bonf)

trans_bonf$FDR = p.adjust(trans_bonf$p.value.bonf, method="fdr")

# head the output
print("Top of trans-by-cis after BH/FDR correction")
trans_bonf = as.data.frame(trans_bonf[order(trans_bonf$FDR),])
head(trans_bonf)

# Save the output
write.table(trans_bonf, sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis_bonf_fdr.txt"), chr, n_expr_pcs))
