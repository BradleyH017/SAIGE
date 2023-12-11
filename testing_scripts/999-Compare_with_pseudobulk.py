#!/usr/bin/env python
# Bradley Dec 2023

# Load packages
import matplotlib as mp
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import scanpy as sc
from pandas_plink import read_plink1_bin
import plotly.express as px
import statsmodels.stats.multitest as smt



# Load the data
catdir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/category__machine/Enterocyte/old2"
outdir=f"{catdir}/plots"
chr="2"
n_expr_pcs="20"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"


# Specify sc/pseudo-bulk files for this chromosome
options = catdir.split("/")
aggregate_on= [element for element in options if "__machine" in element][0]
level = options[-2]
pb=f"/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123/{aggregate_on}-{level}-dMean/cis_nominal1.cis_qtl_pairs.chr{chr}.tsv"
sc=f"{catdir}/chr{chr}_nPC_{n_expr_pcs}.txt"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
geno_dir=f"{general_file_dir}/genotypes"

# Load in the files
pb_res = pd.read_csv(pb, sep = "\t")
sc_res = pd.read_csv(sc, sep = "\t", header=None)
sc_res = sc_res.iloc[:,1:]
sc_res.columns = ["CHR", "POS",  "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA", "Is.SPA", "N", "Gene", "qvalue", "lfdr"]

# Combine the phenotype/variant columns
sc_res['variant_phenotype'] = sc_res['MarkerID'] + "_" + sc_res['Gene']
pb_res['variant_phenotype'] = pb_res['variant_id'] + "_" + pb_res['phenotype_id']

# Find common elements
common = set(pb_res['variant_phenotype']).intersection(sc_res['variant_phenotype'])

# Filter DataFrames for common elements
sc_res = sc_res[sc_res['variant_phenotype'].isin(common)]
pb_res = pb_res[pb_res['variant_phenotype'].isin(common)]

# Order DataFrames by 'variant_phenotype'
sc_res = sc_res.sort_values(by='variant_phenotype')
pb_res = pb_res.sort_values(by='variant_phenotype')

# Perform Pearson correlation test
sc_res['BETA'] = sc_res['BETA'].astype('float')
rho, p_value = pearsonr(pb_res['slope'], sc_res['BETA'])

# Merge DataFrames on 'variant_phenotype'
merged = pd.merge(sc_res, pb_res, on='variant_phenotype')

# Set up Seaborn plot
max_effect = max([np.abs(max(merged['slope'])), np.abs(max(merged['BETA']))])
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
# Scatter plot with hexbin
#sns.jointplot(
#    data=merged, 
#    x='slope', 
#    y='BETA', 
#    kind='hex', 
#    cmap='viridis', 
#    joint_kws=dict(gridsize=70)
#)
sns.scatterplot(x='slope', y='BETA', data=merged)
plt.plot([-max_effect, max_effect], [-max_effect, max_effect], color='red', linestyle='dashed', label='x=y Line')
plt.xlabel('Effect size - TensorQTL')
plt.ylabel('Effect size - SAIGE')
plt.xlim(-max_effect, max_effect)
plt.ylim(-max_effect, max_effect)
plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
plt.title(f"Chromosome {chr} only")
# Save the plot
plt.savefig(f"{outdir}/cor_betas_pb_chr{chr}_py.png")
plt.clf()

# Do the same but now look at p-value
max_effect = max([np.abs(max(-np.log10(merged['pval_nominal']))), np.abs(max(-np.log10(merged['p.value'][merged['p.value'] > 0])))])
merged['tensor_nlog10'] = -np.log10(merged['pval_nominal'])
merged['saige_nlog10'] = -np.log10(merged['p.value'])
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
# Scatter plot with hexbin
#sns.jointplot(
#    data=merged, 
#    x='slope', 
#    y='BETA', 
#    kind='hex', 
#    cmap='viridis', 
#    joint_kws=dict(gridsize=70)
#)
sns.scatterplot(x='tensor_nlog10', y='saige_nlog10', data=merged)
plt.plot([0, max_effect], [0, max_effect], color='red', linestyle='dashed', label='x=y Line')
plt.xlabel('-log10(TensorQTL p-value)')
plt.ylabel('-log10(SAIGEQTL p-value)')
plt.xlim(0, max_effect)
plt.ylim(0, max_effect)
plt.title(f"Chromosome {chr} only")
#plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
# Save the plot
plt.savefig(f"{outdir}/cor_p-values_pb_chr1_py.png")
plt.clf()

# Count the number of times SAIGE and TensorQTL exceed each others significance
all_test = merged.shape[0]
nsaige_win = merged[merged['p.value'] < merged['pval_nominal']].shape[0]
perc_saige_win = "{:.2f}".format(100*(nsaige_win/all_test))
ntensor_win = merged[merged['p.value'] > merged['pval_nominal']].shape[0]
perc_tensor_win = "{:.2f}".format(100*(ntensor_win/all_test))
print(f"Out of a possible {all_test} tests, SAIGE wins {nsaige_win} times ({perc_saige_win}%), Tensor wins {ntensor_win} times ({perc_tensor_win}%)")

# Those with p<0.05 and p<5e-8
cutoffs = [0.05, 5e-8]
for c in cutoffs:
    sub = merged[(merged['p.value'] < c) | (merged['pval_nominal'] < c)]
    all_test = sub.shape[0]
    nsaige_win = sub[sub['p.value'] < sub['pval_nominal']].shape[0]
    perc_saige_win = "{:.2f}".format(100*(nsaige_win/all_test))
    ntensor_win = sub[sub['p.value'] > sub['pval_nominal']].shape[0]
    perc_tensor_win = "{:.2f}".format(100*(ntensor_win/all_test))
    print(f"At a p-value threshold of {c}: Out of a possible {all_test} tests, SAIGE wins {nsaige_win} times ({perc_saige_win}%), Tensor wins {ntensor_win} times ({perc_tensor_win}%)")


# Looking at the top hit per gene (maybe we greater a greater hit rate per gene as opposed to variants?) [Make sure dfs are corrected same way]
pb_res['qvalue_within_gene'] = pb_res.groupby('phenotype_id')['pval_nominal'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
pb_res_min_qvalue_rows = pb_res.loc[pb_res.groupby('phenotype_id')['qvalue_within_gene'].idxmin()]
pb_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(pb_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]
sc_res['qvalue_within_gene'] = sc_res.groupby('Gene')['p.value'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
sc_res_min_qvalue_rows = sc_res.loc[sc_res.groupby('Gene')['qvalue_within_gene'].idxmin()]
sc_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(sc_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]

# Plot the p-values for the most-significant variant / gene
pb_res.rename(columns = {'phenotype_id': 'Gene'}, inplace=True)
merged_top = pd.merge(sc_res_min_qvalue_rows, pb_res_min_qvalue_rows, on='Gene')



# Have some where there is a very large  -log10(p) for TensorQTL, but much smaller for SAIGE. Investigate these
test = merged[(-np.log10(merged.pval_nominal) > 10) & (-np.log10(merged['p.value']) < 10)]
min_index = test.groupby('Gene')['pval_nominal'].idxmin()
test_min_per_gene = test.loc[min_index]
# Save these results
test_min_per_gene.to_csv(f"{outdir}/tensorqtl_beats_saige_top_cases_chr1.csv")

# Load the anndata object to investigate these effects in the cells of interest
adata = ad.read_h5ad(phenotype__file)
adata = adata[adata.obs[aggregate_on] == level]

# Load the genotypes of these results
G = read_plink1_bin(f"{geno_dir}/plink_genotypes_chr{chr}.bed")
G_id = pd.read_csv(f"{geno_dir}/plink_genotypes_chr{chr}.bim", sep = "\t", header=None)
G['variant'] = G_id.iloc[:,1]
G = G[G.sample.isin(adata.obs.Corrected_genotyping_ID)]
G = G[:,G['variant'].isin(test_min_per_gene.MarkerID)]

# Plot the expression of these genes
genotype_df = G.to_dataframe(name='genotype').reset_index()
genotype_df = genotype_df.pivot(index='sample', columns='variant', values='genotype')
expression_df = pd.DataFrame.sparse.from_spmatrix(adata.X)
expression_df.columns=adata.var.index.values
expression_df.index = adata.obs.index.values

# Combine
genotype_df = genotype_df.reset_index()
genotype_df = genotype_df.rename(columns={'sample': 'Corrected_genotyping_ID'})
temp = adata.obs
temp = temp.reset_index()
genotype_df = genotype_df.merge(temp[["index", "Corrected_genotyping_ID"]], on="Corrected_genotyping_ID")
genotype_df = genotype_df.set_index("index")
expression_df = expression_df.merge(genotype_df, left_index=True, right_index=True)

# Make a figure for these results in the single cell expression
for r, marker in enumerate(test_min_per_gene.MarkerID.values):
    gene = test_min_per_gene.Gene.values[r]
    fname=f"{outdir}/{marker}_{gene}.png"
    plt.figure()
    sns.violinplot(x=marker, y=gene, data=expression_df)
    plt.savefig(fname)
    plt.clf()


for r, row in test_min_per_gene.iterrows():
    marker = row['MarkerID']
    print(marker)
    gene = row['Gene']
    fig = px.violin(x=expression_df[marker], y=expression_df[gene], box=True, points="all")
    fig.update_layout(title=f'Violin Plot for MarkerID: {marker}, Gene: {gene}')
    fig.update_xaxes(title_text=marker)
    fig.update_yaxes(title_text=gene)
    fig.update_layout(
        plot_bgcolor='lightgray',  # Set the background color
        paper_bgcolor='white'      # Set the color of the plotting area
        )
    fig.update_traces(marker_color='rgba(100, 100, 255, 0.6)')
    fig.write_image(f"{outdir}/{marker}_{gene}.png")
    

# Really convincing example is: chr1:159965126:G:C on ENSG00000224259 from single cell
# Check the concordance of effects for these very strong cases
min_max = np.abs(np.concatenate((np.array(test_min_per_gene.slope.values) , np.array(test_min_per_gene.BETA.values)))).max()
plt.figure(figsize=(8, 6))
plt.scatter(test_min_per_gene["slope"], test_min_per_gene["BETA"], alpha=0.7)
for i, label in enumerate(test_min_per_gene.variant_phenotype):
    print(label)
    plt.text(test_min_per_gene.slope.values[i], test_min_per_gene.BETA.values[i], label,  ha='center', va='center', rotation=-45)

# Add a red dashed x=y line
plt.plot([-min_max, min_max],
         [-min_max, min_max],
         color='red', linestyle='--', alpha=0.5)
plt.title('Effect sizes of variants with dramatically greater expression in TensorQTL vs SAIGE')
plt.xlabel('TensorQTL beta (large significance variants)')
plt.ylabel('SAIGEQTL beta (large effect variants)')
plt.savefig(f"{outdir}/cor_beta_large_effect_variants.png", bbox_inches='tight')
plt.clf()

