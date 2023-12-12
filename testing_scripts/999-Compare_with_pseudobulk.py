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
from scipy.stats import linregress
from scipy import stats
import scipy.stats as stats
import statsmodels.api as sm 
from matplotlib_venn import venn2





# Load the data (multiple chromosomes)
catdir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/label__machine/T_cell_CD8_1"
outdir=f"{catdir}/plots"
chr=range(1,2)
formatted_range = f"{chr[0]}_{chr[-1]}"
n_expr_pcs="15"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"


# Specify sc/pseudo-bulk files for this chromosome
options = catdir.split("/")
aggregate_on= [element for element in options if "__machine" in element][0]
level = options[-1]
pb_base=f"/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123/{aggregate_on}-{level}-dMean"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
geno_dir=f"{general_file_dir}/genotypes"

# Load in the pseudobulk
pb_res = []
for c in chr:
    print(c)
    fname = f"{pb_base}/cis_nominal1.cis_qtl_pairs.chr{str(c)}.tsv"
    f = pd.read_csv(fname, sep = "\t")
    pb_res.append(f)

pb_res = pd.concat(pb_res, axis=0, ignore_index=True)

# Load in the sc_res
sc_res = []
for c in chr:
    print(c)
    fname = f"{catdir}/chr{str(c)}_nPC_{n_expr_pcs}.txt"
    f = pd.read_csv(fname, sep = "\t", header=None, usecols=range(17))
    sc_res.append(f)

sc_res = pd.concat(sc_res, axis=0, ignore_index=True)
sc_res.columns = ["CHR", "POS",  "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA", "Is.SPA", "N", "Gene"]
sc_res = sc_res.dropna(subset=['Gene'])
mask = sc_res['Gene'].notna()
sc_res = sc_res[mask & sc_res['Gene'].str.contains("ENSG")]
sc_res['MarkerID'] = sc_res['MarkerID'].astype(str).replace('nan', '')
sc_res['Gene'] = sc_res['Gene'].astype(str).replace('nan', '')

# Combine the phenotype/variant columns
sc_res['variant_phenotype'] = sc_res['MarkerID'] + "_" + sc_res['Gene']
sc_res = sc_res.drop_duplicates(subset=['variant_phenotype'])
pb_res['variant_phenotype'] = pb_res['variant_id'] + "_" + pb_res['phenotype_id']

# Find common tests
common = set(pb_res['variant_phenotype']).intersection(sc_res['variant_phenotype'])

# Before subset - plot histogram of p-values in SAIGE
sc_res['p.value'] = sc_res['p.value'].astype('float')
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
plt.hist(sc_res['p.value'], bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
# Add labels and title
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of SAIGE p-alues')
plt.savefig(f"{outdir}/SAIGE_pvalues_chr{formatted_range}_py.png")
plt.clf()

# Also QQ plot
plt.figure(figsize=(8, 6))
obs = sc_res['p.value']
n = len(obs)
obs_p = -np.log10(np.sort(obs))
th_p = np.arange(1/float(n), 1 + (1/float(n)), 1/float(n))
th_p = -np.log10(th_p)
# Scatter plot with small blue points
plt.scatter(th_p, obs_p, color='blue', s=2)
# Red dashed line
x = np.linspace(*plt.xlim())
plt.plot(x, x, color='red', linestyle='--')
# Adjustments
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.title("QQ plot - Observed vs Expected p-value Distribution (-log10) - All Tests")
plt.tight_layout()  # Removes white margin between the plot and axis
# Save or show the plot
plt.savefig(f"{outdir}/SAIGE_pvalues_chr{formatted_range}_qq_py.png")

# Do the same but for three random genes
random_genes = np.random.choice(sc_res['Gene'].unique(), 3, replace=False)
plt.figure(figsize=(8, 6))
for gene in random_genes:
    # Filter data for the current gene
    gene_data = sc_res[sc_res['Gene'] == gene]['p.value']
    n = len(gene_data)
    obs_p = -np.log10(np.sort(gene_data))
    th_p = -np.log10(np.arange(1/float(n), 1 + (1/float(n)), 1/float(n)))
    # Scatter plot with small points and different colors
    plt.scatter(th_p, obs_p, label=f'Gene: {gene}', s=2)

# Red dashed line
x = np.linspace(*plt.xlim())
plt.plot(x, x, color='red', linestyle='--')
# Adjustments
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.title("QQ plot - Observed vs Expected p-value Distribution (-log10) - Random Genes")
plt.legend()  # Show legend with gene labels
plt.tight_layout()
# Save or show the plot
plt.savefig(f"{outdir}/SAIGE_pvalues_chr{formatted_range}_qq_py_random_genes.png")



# Filter DataFrames for common elements
sc_res = sc_res[sc_res['variant_phenotype'].isin(common)]
pb_res = pb_res[pb_res['variant_phenotype'].isin(common)]

# Order DataFrames by 'variant_phenotype'
sc_res = sc_res.sort_values(by='variant_phenotype')
pb_res = pb_res.sort_values(by='variant_phenotype')

# Perform Pearson correlation test and linear regression
sc_res['BETA'] = sc_res['BETA'].astype('float')
rho, p_value = pearsonr(pb_res['slope'], sc_res['BETA'])
slope, intercept, r_value, p_value, std_err = linregress(pb_res['slope'], sc_res['BETA'])


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
#x_line = np.linspace(pb_res['slope'], pb_res['slope'], 100)
#y_line = slope * x_line + intercept
# plt.plot(x_line, y_line, color='red', label=f'Linear Model: y = {slope:.2f}x + {intercept:.2f}')
plt.plot([-max_effect, max_effect], [-max_effect, max_effect], color='red', linestyle='dashed', label='x=y Line')
plt.xlabel('Effect size - TensorQTL')
plt.ylabel('Effect size - SAIGE')
plt.xlim(-max_effect, max_effect)
plt.ylim(-max_effect, max_effect)
plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
plt.title(f"Chromosome {formatted_range} only")
# Save the plot
plt.savefig(f"{outdir}/cor_betas_pb_chr{formatted_range}_py.png")
plt.clf()

# Do the same but now look at p-value
max_effect = max([np.abs(max(-np.log10(merged['pval_nominal']))), np.abs(max(-np.log10(merged['p.value'].astype('float'))))])
merged['tensor_nlog10'] = -np.log10(merged['pval_nominal'])
merged['saige_nlog10'] = -np.log10(merged['p.value'].astype(float))
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
plt.title(f"Chromosome {formatted_range} only")
#plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
# Save the plot
plt.savefig(f"{outdir}/cor_p_values_pb_chr{formatted_range}_py.png")
plt.clf()

# Count the number of times SAIGE and TensorQTL exceed each others significance
merged['p.value'] = merged['p.value'].astype('float')
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

# Plotting the overlap of genes detected by each method at different resolutions
pb_res_min_qvalue_rows.rename(columns={"phenotype_id": "Gene"}, inplace=True)
merged_qvalue = sc_res_min_qvalue_rows.merge(pb_res_min_qvalue_rows, on="Gene", suffixes=('_saige', '_tensor')) 
cutoffs = [1, 0.05]
for c in cutoffs:
    sub = merged_qvalue[(merged_qvalue['qvalue_within_gene_saige'] < c) | (merged_qvalue['qvalue_within_gene_tensor'] < c)]
    all_test = sub.shape[0]
    nsaige_win = sub[sub['qvalue_within_gene_saige'] < sub['qvalue_within_gene_tensor']].shape[0]
    perc_saige_win = "{:.2f}".format(100*(nsaige_win/all_test))
    ntensor_win = sub[sub['qvalue_within_gene_saige'] > sub['qvalue_within_gene_tensor']].shape[0]
    perc_tensor_win = "{:.2f}".format(100*(ntensor_win/all_test))
    print(f"At a p-value threshold of {c}: Out of a possible {all_test} tests, SAIGE wins {nsaige_win} times ({perc_saige_win}%), Tensor wins {ntensor_win} times ({perc_tensor_win}%)")

# Plot the overlap of most interesting results
sub = merged_qvalue[(merged_qvalue['qvalue_within_gene_saige'] < 0.05) | (merged_qvalue['qvalue_within_gene_tensor'] < 0.05)]
saige_hits = set(sub[sub['qvalue_within_gene_saige'] < 0.05]['Gene'].tolist())
tensor_hits = set(sub[sub['qvalue_within_gene_tensor'] < 0.05]['Gene'].tolist())
venn_labels = {'100': 'SAIGE-QTL', '010': 'TensorQTL', '110': 'Intersection'}

plt.figure(figsize=(8, 6))
venn2(subsets=[saige_hits, tensor_hits], set_labels=('SAIGE-QTL', 'TensorQTL'), set_colors=('skyblue', 'lightgreen'), alpha=0.7)
plt.savefig(f"{outdir}/venn_saige_tensor_fdr_0.05_chr{formatted_range}_py.png")
plt.clf()

# Also plot the correlation of betas/p-values for these significant effects - if only significant in one condition, plot one colour and another if the opposite
max_effect = max([np.abs(max(sub['slope'])), np.abs(max(sub['BETA']))])
sub['exclusivity'] = 'both'  # Default value for all rows
sub.loc[(sub['qvalue_within_gene_saige'] < 0.05) & (sub['qvalue_within_gene_tensor'] > 0.05), 'exclusivity'] = 'SAIGE_only'
sub.loc[(sub['qvalue_within_gene_saige'] > 0.05) & (sub['qvalue_within_gene_tensor'] < 0.05), 'exclusivity'] = 'Tensor_only'
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='slope', y='BETA', hue='exclusivity', data=sub)
#x_line = np.linspace(pb_res['slope'], pb_res['slope'], 100)
#y_line = slope * x_line + intercept
# plt.plot(x_line, y_line, color='red', label=f'Linear Model: y = {slope:.2f}x + {intercept:.2f}')
plt.plot([-max_effect, max_effect], [-max_effect, max_effect], color='red', linestyle='dashed', label='x=y Line')
plt.xlabel('Effect size - TensorQTL')
plt.ylabel('Effect size - SAIGE')
plt.xlim(-max_effect, max_effect)
plt.ylim(-max_effect, max_effect)
plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
plt.title(f"FDR < 0.05, Chromosome {formatted_range} only")
# Save the plot
plt.savefig(f"{outdir}/cor_betas_pb_fdr_0.05_one_both_chr{formatted_range}_py.png")
plt.clf()

# Same but p-values
max_effect = max([np.abs(max(-np.log10(sub['pval_nominal']))), np.abs(max(-np.log10(sub['p.value'].astype('float'))))])
sub['tensor_nlog10'] = -np.log10(sub['pval_nominal'])
sub['saige_nlog10'] = -np.log10(sub['p.value'].astype(float))
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='tensor_nlog10', y='saige_nlog10', data=sub, hue='exclusivity')
plt.plot([0, max_effect], [0, max_effect], color='red', linestyle='dashed', label='x=y Line')
plt.xlabel('-log10(TensorQTL p-value)')
plt.ylabel('-log10(SAIGEQTL p-value)')
plt.xlim(0, max_effect)
plt.ylim(0, max_effect)
plt.title(f"FDR < 0.05, Chromosome {formatted_range} only")
#plt.text(5.5, 0.5, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
# Save the plot
plt.savefig(f"{outdir}/cor_p_values_pb_fdr_0.05_chr{formatted_range}_py.png")
plt.clf()

# For nominal p < 5e-8
sub = merged_qvalue[(merged_qvalue['p.value'] < 5e-8) | (merged_qvalue['pval_nominal'] < 5e-8)]
saige_hits = set(sub[sub['qvalue_within_gene_saige'] < 5e-8]['Gene'].tolist())
tensor_hits = set(sub[sub['qvalue_within_gene_tensor'] < 5e-8]['Gene'].tolist())
venn_labels = {'100': 'SAIGE-QTL', '010': 'TensorQTL', '110': 'Intersection'}
plt.figure(figsize=(8, 6))
venn2(subsets=[saige_hits, tensor_hits], set_labels=('SAIGE-QTL', 'TensorQTL'), set_colors=('skyblue', 'lightgreen'), alpha=0.7)
plt.savefig(f"{outdir}/venn_saige_tensor_p_5e-8_chr{formatted_range}_py.png")
plt.clf()


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

