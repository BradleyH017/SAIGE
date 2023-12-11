#!/usr/bin/env Rscript
# Bradley Nov 2023

# Inherit options
option_list = list(
    make_option(c("-c", "--catdir"), action = "store", default = NA, type ="character",
        help="Where are the files located?"),
    make_option(c("-chr", "--chromosome"), action = "store", default = NA, type ="character",
        help="Which chromosome to test"),
    make_option(c("-n", "--n_expr_pcs"), action = "store", default = NA, type ="character",
        help="What is the optimum number of expression PCs to use")
    )
opt = parse_args(OptionParser(option_list = option_list))
catdir = opt$c
chr= opt$chr
n_expr_pcs = opt$n

# testing
catdir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/category__machine/Enterocyte/old"
outdir=paste0(catdir, "/plots")

# Load libraries
library(Cairo)
library(ggplot2)

# Specify sc/pseudo-bulk files for this chromosome
options = unlist(strsplit(catdir, "\\/"))
aggregate_on=options[grep("machine", options)]
level = options[length(options)-1]
pb=paste0("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123/", aggregate_on, "-", level, "-dMean/cis_nominal1.cis_qtl_pairs.chr", chr, ".tsv")
sc=paste0(catdir, "/chr", chr, "_nPC_", n_expr_pcs, ".txt", sep = "")

# Load in each (pb is not gzipped, sc is)
pb_res = read.delim(pb, header=T)
sc_res = read.delim(sc, header=F)

# Add the header for the sc results
sc_res_header=c("CHR", "POS",  "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA", "Is.SPA", "N", "Gene")
colnames(sc_res) = sc_res_header

# Combine the phenotype/variant columns on each file and filter for those in common, then order by these
sc_res$variant_phenotype = paste0(sc_res$MarkerID, "_", sc_res$Gene)
pb_res$variant_phenotype = paste0(pb_res$variant_id, "_", pb_res$phenotype_id)
common = intersect(pb_res$variant_phenotype, sc_res$variant_phenotype)
sc_res = sc_res[sc_res$variant_phenotype %in% common,]
pb_res = pb_res[pb_res$variant_phenotype %in% common,]
sc_res = sc_res[order(sc_res$variant_phenotype),]
pb_res = pb_res[order(pb_res$variant_phenotype),]

# Now look at the correlation between the betas
b_cor = cor.test(pb_res$slope, as.numeric(sc_res$BETA))
if(file.exists(outdir) != T){
    dir.create(outdir)
}
max_effect_size = max(abs(pb_res$slope), abs(as.numeric(sc_res$BETA)))
sc_res$BETA = as.numeric(sc_res$BETA)
label=paste0("Pearsons r=", signif(b_cor$estimate,2), "\np=", signif(b_cor$p.value, 2))
merged = merge(sc_res, pb_res, by="variant_phenotype")
beta_density_plot <- ggplot(merged, aes(x = slope, y = BETA)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() + 
  labs(title = "Concordance of effect sizes (chr1)", x = "TensorQTL", y = "SAIGEQTL") +
  # Add label in the top-left corner
  annotate("text", x = -max_effect_size, y = max_effect_size, label = label, vjust = 1, hjust = 0, color = "black", size = 4) +
  # Add red dashed line for x=y
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  geom_smooth(method = "lm", se = FALSE, color = "blue", formula = y ~ x, linetype = "solid") + 
  annotate("text", x = 0.5, y = -0.9, label = paste("Equation: y = ", round(coef(lm(BETA ~ slope, data = merged))[2], 2), "x + ", round(coef(lm(BETA ~ slope, data = merged))[1], 2), "\nR-squared: ", round(summary(lm(BETA ~ slope, data = merged))$r.squared, 3)), color = "blue", size = 4)

ggsave(paste0(outdir, "/cor_betas_pb_chr1.png"), plot = beta_density_plot, width = 8, height = 6, units = "in")


p_plot=ggplot(merged, aes(-log10(pval_nominal), -log10(p.value))) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Concordance of p-val sizes (chr1)", x = "-log10(TensorQTL p)", y = "-log10(SAIGEQTL p") +
  # Add red dashed line for x=y
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # Customize theme (optional)
  theme_bw()

ggsave(paste0(outdir, "/cor_pval_pb_chr1.png"), plot = p_plot, width = 8, height = 6, units = "in")

# Can see some instances where the p-values for TensorQTL far surpass those of SAIGE. 
# Investigating these
test = merged[-log10(merged$pval_nominal) > 10 & -log10(merged$p.value) < 10, ]