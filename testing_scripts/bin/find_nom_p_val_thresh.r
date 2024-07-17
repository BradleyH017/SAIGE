#!/usr/bin/env Rscript

# Bradley Jan 2024
library(optparse)

# Inherit options
option_list = list(
    make_option(c("-c", "--chrom_wide_file"), action = "store", default = NA, type ="character",
        help="chromosome_wide_file"), 
    make_option(c("-f", "--file"), action = "store", default = NA, type ="character",
        help="File to find p-val thresh for"), 
    make_option(c("-t", "--threshold"), action = "store", default = NA, type ="character",
        help="Q value threshold")
    )
opt = parse_args(OptionParser(option_list = option_list))
chrom_file=opt$c
file = opt$f
threshold = as.numeric(opt$t)

# Load in the chromosome wide file
chrom = read.delim(chrom_file)

# Find minimum qvalue for gene (corrected across genes)
fname = unlist(strsplit(file, "\\/"))
fname = fname[length(fname)]
gene = unlist(strsplit(fname, "\\_"))[1]
chrom = chrom[!is.na(chrom[,17]),]
min_q = chrom[chrom[,17] == gene,]$qvalues_across_genes
if(min_q > threshold){
    stop("No cis effect detected across chromosoomes, stopping")
}

# Load in actual file
res = read.delim(file)

# Subset sig (below this across gene threshold)
res_sig = res[as.numeric(res$qvalues) < min_q,]

# Otherwise, extract maximum p-value at this threshold
max_pval = max(res_sig$`p.value`)

# Save this
fname = unlist(strsplit(file, "\\/"))
fname = fname[length(fname)]
basedir = gsub(fname, "", file)
cond_dir = paste0(basedir, "conditional")
base_fname = gsub(".txt", "", fname)
write.table(max_pval, paste0(cond_dir, "/", base_fname, "_pvalue_thresh.txt"), quote=F, row.names=F, col.names=F)