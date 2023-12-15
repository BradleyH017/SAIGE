#!/usr/bin/env Rscript

# Bradley Dec 2023

library(optparse)
library(dplyr)
library(qvalue)

# Inherit options
option_list = list(
    make_option(c("-f", "--file"), action = "store", default = NA, type ="character",
        help="File to estimate q-values for"),
    make_option(c("-c", "--column"), action = "store", default = NA, type ="character",
        help="Which column to calculate p-value correction of"),
    make_option(c("-n", "--new_q_val_column"), action = "store", default = NA, type ="character",
        help="What do you want to call the new column"),
    make_option(c("-w", "--within_gene"), action = "store", default = NA, type ="character",
        help="within genes or across genes?")
    )
opt = parse_args(OptionParser(option_list = option_list))
file = opt$f
column = as.numeric(opt$c)
new_column = opt$n
within = opt$w

# Load data
res = read.delim(file, header=F)
if(is.na(as.numeric(res[,as.numeric(column)][1]))){
    colnames(res) = res[1,]
    res = res[-1,]
}

# Extract p_values
p_values = as.numeric(res[,as.numeric(column)])

# Generate q-value object
qobj = qvalue(p = p_values)

# Add q-values
res[,new_column] = qobj$qvalues

# Add local fdr
res$lfdr = qobj$lfdr

# Replace the original file with the current one
write.table(res, file, col.names=T, quote=F, sep = "\t",row.names=F)

# Save a file with the variant with the minimum qvalue for this gene
min_file = gsub(".txt", "", file)
min_file = paste0(min_file, "_minimum_q.txt")
to_save = res[res[,new_column] == min(res[,new_column]),]
# If little significance, may get lots of results with the minimum qvalue
if(nrow(to_save > 1)){
    to_save=res[res[,new_column] == min(res[,new_column]),]
    # If have variants in complete LD - take the first
    if(nrow(to_save) > 1){
        to_save = to_save[1,]
    }
}

if(within == "TRUE"){
    write.table(to_save, min_file, col.names=T, quote=F, sep = "\t", row.names=F)
}
print("Performed q-value correction")