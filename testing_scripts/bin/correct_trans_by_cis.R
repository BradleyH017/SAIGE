#!/usr/bin/env Rscript

###### Bradley Dec 2023
library(optparse)
library(dplyr)
library(ggplot2)

# Inherit options
args <- commandArgs(trailingOnly = TRUE)
cat("Command line arguments:", paste(args, collapse = " "), "\n")
# Parse command-line arguments
catdir_index <- which(args == "-d")
chr_index <- which(args == "-c")
n_expr_pcs_index <- which(args == "-n")
window_index <- which(args == "-w")
annotation_file_index <- which(args == "-a")
catdir <- args[catdir_index + 1]
chr <- args[chr_index + 1]
n_expr_pcs <- args[n_expr_pcs_index + 1]
window <- args[window_index + 1]
annotation_file <- args[annotation_file_index + 1]
# Print them
cat("catdir:", catdir, "\n")
cat("chr:", chr, "\n")
cat("n_expr_pcs:", n_expr_pcs, "\n")
cat("window:", window, "\n")
cat("annotation_file:", annotation_file, "\n")
print(sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis.txt"), as.character(chr), as.character(n_expr_pcs)))

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
print("Filtering trans to remove out cis")
filtered_trans =
  trans %>%
  filter(CHR != gene_chromosome | 
         (gene_strand == '-' & abs(as.numeric(POS) - as.numeric(gene_end)) >= as.numeric(window)) |
         (gene_strand == '+' & abs(as.numeric(POS) - as.numeric(gene_start)) >= as.numeric(window)))

# Perform BH correction within genes
print("Performing BH correction within gene")
filtered_trans <- filtered_trans %>%
  group_by(Gene) %>%
  mutate(p.value.bonf = p.adjust(p.value, method = "bonferroni")) %>% as.data.frame()

# Plot a distribution of the nominal p-values 
histogram <- ggplot(filtered_trans, aes(x = p.value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Nominal p-values (trans-by-cis)", x = "Values", y = "Frequency")

# Save the plot using ggsave
#ggsave(sprintf("../plots/histogram_trans_p.png"), plot = histogram, width = 8, height = 6, dpi = 300)

# Save the filtered output with bonferoni
write.table(filtered_trans, sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis.txt"), as.character(chr), as.character(n_expr_pcs)), col.names=T, quote=F, sep = "\t", row.names=F)

# Subset for those effects that are minimum per gene
filtered_trans <- filtered_trans %>%
  group_by(Gene) %>%
  filter(rank(p.value.bonf, ties.method = "min") == 1) %>%
  arrange(p.value) %>%
  slice(1) %>%
  ungroup() %>% as.data.frame()

# head the output
print("Top of trans-by-cis after BH/FDR correction")
filtered_trans = as.data.frame(filtered_trans[order(filtered_trans$`p.value.bonf`),])
head(filtered_trans)

# Plot histogram of bonferoni p-values for this chromosome
histogram <- ggplot(filtered_trans, aes(x = p.value.bonf)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Bonferoni p-values (trans-by-cis)", x = "Values", y = "Frequency")
#ggsave(sprintf("../plots/histogram_trans_p_bonf.png"), plot = histogram, width = 8, height = 6, dpi = 300)

# Save the output
write.table(filtered_trans, sprintf(paste0(catdir, "/chr%s_nPC_%s_trans_by_cis_no_cis_bonf.txt"), chr, n_expr_pcs), col.names=T, quote=F, sep = "\t", row.names=F)
