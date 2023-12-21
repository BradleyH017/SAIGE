#!/usr/bin/env Rscript

###### Bradley Nov 2023
library(optparse)
library(dplyr)


# Inherit options
option_list = list(
    make_option(c("-c", "--catdir"), action = "store", default = NA, type ="character",
        help="Where are the files located?")
    )
opt = parse_args(OptionParser(option_list = option_list))
catdir = opt$c

# Load files and process
group_files=list.files(catdir)
group_files=group_files[grep("_minimum_q_all_genes", group_files)]
nPCs = gsub("_minimum_q_all_genes.txt", "", group_files)
nPCs = gsub("chr1_nPC_", "", nPCs)

res = data.frame("nExpr_PCs"=nPCs, "neGenes"="")
for(f in nPCs){
    temp = read.delim(paste0(catdir, "/chr1_nPC_", f, "_minimum_q_all_genes.txt"), header=T)
    nsig = nrow(temp[temp$qvalues_across_genes < 0.05,])
    res$neGenes[which(nPCs == f)] = nsig
    rm(temp)
}
res=res[order(as.numeric(res$nExpr_PCs)),]
print("Number of eGenes per conditions (FDR < 0.05)")
print(res)
write.table(res, paste0(catdir, "/expr_nPCs_check_chr1.txt"), row.names=F, col.names=F, quote=F)

# Save the value of the nPCs at which the maximum number of 
max_res = res[res$neGenes == max(as.numeric(res$neGenes)),]$nExpr_PCs
if(length(max_res) > 1){
    max_res = max_res[1]
}
write.table(max_res, paste0(catdir, "/optim_nPCs_chr1.txt"), row.names=F, col.names=F, quote=F)
