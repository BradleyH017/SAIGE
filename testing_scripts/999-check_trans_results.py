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
from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import plotly.graph_objects as go
from matplotlib.lines import Line2D



# Load in the trans results
label="Pericytes"
n_expr_pcs=5
catdir = f"/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/label__machine/{label}"
res = pd.read_csv(f"{catdir}/trans/all_nPC_{str(n_expr_pcs)}_trans_by_cis_no_cis_bonf.txt", sep = "\t")
geno_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/genotypes"
outdir=f"{catdir}/plots"

# Grab the top 5 and plot their trans-effects
res_sig = res[res['bonf_qvalue_across_genes'] < 0.01]
res_sig = res_sig.sort_values(by='p.value.bonf')
res_plot = res_sig.head(20)

# Plot distribution of significant betas
sns.distplot(res_sig['BETA'], hist=True, rug=True)
plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.title('Distribution effect size estimates')
plt.xlabel('Effect size')
#plt.axvline(x = 0.7, color = 'red', linestyle = '--', alpha = 0.5)
plt.savefig(f"{outdir}/sig_trans_effect_sizes.png", bbox_inches='tight')
plt.clf()

# Make a circos plot of the significant effects (https://moshi4.github.io/pyCirclize/circos_plot/)
chr_bed_file, cytoband_file, chr_links = load_eukaryote_example_dataset("hg38")
# Initialize Circos from BED chromosomes
circos = Circos.initialize_from_bed(chr_bed_file, space=3)
circos.text(f"{label}\ntrans-eQTLs", deg=315, r=150, size=12)
circos.add_cytoband_tracks((95, 100), cytoband_file)
# Create chromosome color dict
ColorCycler.set_cmap("hsv")
chr_names = [s.name for s in circos.sectors]
colors = ColorCycler.get_color_list(len(chr_names))
chr_name2color = {name: color for name, color in zip(chr_names, colors)}

# Plot chromosome name & xticks
for sector in circos.sectors:
    sector.text(sector.name, r=120, size=10, color=chr_name2color[sector.name])
    sector.get_track("cytoband").xticks_by_interval(
        40000000,
        label_size=8,
        label_orientation="vertical",
        label_formatter=lambda v: f"{v / 1000000:.0f} Mb",
    )

# Redefine the chr_links using our own data
temp = res_sig[['CHR', 'POS', 'gene_chromosome', 'gene_start', 'gene_end', 'BETA']]
temp.columns = ['ref_chr', 'ref_start', 'query_chr', 'query_start', 'query_end', 'BETA']
temp['ref_chr'] = temp['ref_chr'].apply(lambda x: 'chr' + str(x))
temp['query_chr'] = temp['query_chr'].apply(lambda x: 'chr' + str(x))
temp['ref_end'] = temp['ref_start'].astype('float') + 1
# Determine the color scale based on 'BETA'
min_beta = temp['BETA'].min()
max_beta = temp['BETA'].max()
# Define a colormap based on the minimum and maximum BETA values
cmap = plt.cm.coolwarm  # You can choose other colormaps
norm = plt.Normalize(min_beta, max_beta)
chr_name2color = {name: cmap(norm(value)) for name, value in zip(chr_names, temp['BETA'].astype('float'))}

temp.sort_index(axis=1, inplace=True)
class ChrLink:
    def __init__(self, query_chr, query_start, query_end, ref_chr, ref_start, ref_end, BETA):
        self.query_chr = query_chr
        self.query_start = query_start
        self.query_end = query_end
        self.ref_chr = ref_chr
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.beta = BETA
        
    def __repr__(self):
        return f"ChrLink(query_chr='{self.query_chr}', query_start={self.query_start}, query_end={self.query_end}, ref_chr='{self.ref_chr}', ref_start={self.ref_start}, ref_end={self.ref_end}, beta={self.beta})"

# Convert DataFrame rows to ChrLink objects
chrlink_objects = [ChrLink(**row) for index, row in temp.iterrows()]

# Make colour scale
cmap = plt.cm.coolwarm  # You can choose other colormaps
norm = TwoSlopeNorm(vmin=float(min_beta), vcenter=0, vmax=float(max_beta))
chrlink_colors = [cmap(norm(value)) for value in temp['BETA'].astype('float')]

# Plot chromosome link
for i, link in enumerate(chrlink_objects):
    region1 = (link.query_chr, float(link.query_start), float(link.query_end))
    region2 = (link.ref_chr, float(link.ref_start), float(link.ref_end))
    color = chrlink_colors[i]
    circos.link(region1, region2, color=color, direction=1, alpha=1)

fig = circos.plotfig()
fig.savefig(f"{outdir}/trans_circos_plot.png")

# Load the anndata object
adata = ad.read_h5ad("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad")
options = catdir.split("/")
aggregate_on= [element for element in options if "__machine" in element][0]
level = options[-1]
adata = adata[adata.obs[aggregate_on] == level]

# Load the genotypes of these results
G = read_plink1_bin(f"{geno_dir}/plink_genotypes_cis_{level}.bed")
G_id = pd.read_csv(f"{geno_dir}/plink_genotypes_cis_{level}.bim", sep = "\t", header=None)
G['variant'] = G_id.iloc[:,1]
G = G[G.sample.isin(adata.obs.Corrected_genotyping_ID)]
G = G[:,G['variant'].isin(res_plot.MarkerID)]

# Plot the expression of these genes
genotype_df = G.to_dataframe(name='genotype').reset_index()
genotype_df = genotype_df.pivot(index='sample', columns='variant', values='genotype')
expression_df = pd.DataFrame.sparse.from_spmatrix(adata.X)
expression_df.columns=adata.var.index.values
expression_df.index = adata.obs.index.values
expression_df = expression_df.iloc[:,expression_df.columns.isin(res_plot.Gene)]

# Combine
genotype_df = genotype_df.reset_index()
genotype_df = genotype_df.rename(columns={'sample': 'Corrected_genotyping_ID'})
temp = adata.obs
temp = temp.reset_index()
genotype_df = genotype_df.merge(temp[["index", "Corrected_genotyping_ID"]], on="Corrected_genotyping_ID")
genotype_df = genotype_df.set_index("index")
expression_df = expression_df.merge(genotype_df, left_index=True, right_index=True)
expression_df['sample'] = expression_df.index.str.split('-').str[2]
expression_df['sample'] = expression_df['sample'].astype(str)

# Make a figure for these results in the single cell expression.
# Results are always given for the minor allele with SAIGE
# The homozygous for the minor allele in the genotypes is encoded as '0'
def replace_with_ref_nonref(value, ref, nonref):
    if value == '0.0':
        return nonref * 2  # Repeat 'nonref' twice for '0'
    elif value == '1.0':
        return ref + nonref  # Paste 'ref' and 'nonref' for '1'
    elif value == '2.0':
        return ref * 2  # Repeat 'ref' twice for '2'
    return value

for r, row in res_plot.iterrows():
    marker = row['MarkerID']
    gene = row['Gene']
    ref = row['Allele1']
    nonref = row['Allele2']
    # Create a new DataFrame to store the data for the violin plot
    violin_df = expression_df.copy()
    # Filter rows where 'marker' column is not NaN
    #violin_df = violin_df[violin_df[marker].notna()]
    # Recode the variants (0=minor=nonref)
    violin_df[marker] = violin_df[marker].astype('str')
    values_series = pd.Series(violin_df[marker].values.flatten())
    replace_vectorized = np.vectorize(lambda x: replace_with_ref_nonref(x, ref, nonref))
    violin_df['use_geno'] = replace_vectorized(values_series)
    violin_df = violin_df[violin_df['use_geno'].notna()]
    # Calculate marker counts
    marker_counts = violin_df['use_geno'].value_counts().to_dict()
    # Modify the 'X_Axis' column to include the desired labels
    violin_df['X_Axis'] = (
        violin_df['use_geno'].astype(str) + ' (' +
        violin_df['use_geno'].map(lambda x: str(marker_counts[x])) + ' cells, ' +
        violin_df.groupby('use_geno')['sample'].transform('nunique').astype(str) + ' samples)'
    )
    # Remove the nan
    violin_df = violin_df[~violin_df['X_Axis'].str.contains('nan')]
    # Sort 'X_Axis' labels with custom sorting function
    plot_levels = np.unique(violin_df['X_Axis'])
    substring_ref_2 = ref * 2
    substring_ref_nonref = ref + nonref
    substring_nonref_2 = nonref * 2
    # Filter the plot_levels array for each category
    plot_levels_ref_2 = [level for level in plot_levels if substring_ref_2 in level]
    plot_levels_ref_nonref = [level for level in plot_levels if substring_ref_nonref in level]
    plot_levels_nonref_2 = [level for level in plot_levels if substring_nonref_2 in level]
    plot_levels_ordered = plot_levels_ref_2 + plot_levels_ref_nonref + plot_levels_nonref_2
    violin_df['X_Axis'] = pd.Categorical(violin_df['X_Axis'], categories=plot_levels_ordered, ordered=True)
    # Create the violin plot with custom x-axis labels
    fig = px.violin(
        x=violin_df['X_Axis'],
        y=violin_df[gene],
        box=True,
        points="all"
    )
    # Customize the plot layout
    fig.update_layout(title=f'Violin Plot for MarkerID: {marker}, Gene: {gene}')
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_xaxes(title_text=f'{marker} (Count)')
    fig.update_yaxes(title_text=gene)
    fig.update_layout(
        plot_bgcolor='white',  # Set the background color
        paper_bgcolor='white'  # Set the color of the plotting area
    )
    fig.update_traces(marker_color='rgba(100, 100, 255, 0.6)')
    # Draw a line and annotate with 'BETA' value from the current row
    x_values = np.array([int(value.split()[0]) for value in violin_df[marker].astype('float')])
    y_values = np.array([violin_df[violin_df['X_Axis'] == level][gene].mean() for level in np.unique(violin_df[marker])])
    slope, intercept = np.polyfit(x_values, y_values, 1)
    # Plot linear regression line
    fig.add_trace(go.Scatter(
        x=x_values,
        y=slope * x_values + intercept,
        mode='lines',
        name='Linear Regression',
        line=dict(color='red', width=2)
    ))
    # Add text with 'BETA' and 'bonf_qvalue_across_genes' values
    fig.add_annotation(
        text=f'Linear Regression: y = {slope:.2f}x + {intercept:.2f}',
        x=0.5,
        y=1.1,  # Adjust the y-coordinate to move the text to the top
        showarrow=False,
        font=dict(size=10),
        xanchor='center',  # Center the text horizontally
        yanchor='bottom'   # Align the text to the top
    )
    fig.update_xaxes(tickmode='auto')
    fig.update_yaxes(tickmode='auto')
    # Save the plot to an image file
    fig.write_image(f"{outdir}/{marker}_{gene}_trans.png")

for r, row in res_plot.iterrows():
    try:
        marker = row['MarkerID']
        gene = row['Gene']
        ref = row['Allele1']
        nonref = row['Allele2']
        # Create a new DataFrame to store the data for the violin plot
        violin_df = expression_df.copy()
        order = ['2.0', '1.0', '0.0'] # order x axis
        violin_df[marker] = violin_df[marker].astype('str')
        values_series = pd.Series(violin_df[marker].values.flatten())
        replace_vectorized = np.vectorize(lambda x: replace_with_ref_nonref(x, ref, nonref))
        violin_df['use_geno'] = replace_vectorized(values_series)
        violin_df = violin_df[violin_df['use_geno'].notna()]
        # Calculate marker counts
        marker_counts = violin_df['use_geno'].value_counts().to_dict()
        # Modify the 'X_Axis' column to include the desired labels
        violin_df['X_Axis'] = (
            violin_df['use_geno'].astype(str) + '\n(' +
            violin_df['use_geno'].map(lambda x: str(marker_counts[x])) + ' cells,\n' +
            violin_df.groupby('use_geno')['sample'].transform('nunique').astype(str) + ' samples)'
        )
        # Remove the nan
        violin_df = violin_df[~violin_df['X_Axis'].str.contains('nan')]
        # Sort 'X_Axis' labels with custom sorting function
        plot_levels = np.unique(violin_df['X_Axis'])
        substring_ref_2 = ref * 2
        substring_ref_nonref = ref + nonref
        substring_nonref_2 = nonref * 2
        # Filter the plot_levels array for each category
        plot_levels_ref_2 = [level for level in plot_levels if substring_ref_2 in level]
        plot_levels_ref_nonref = [level for level in plot_levels if substring_ref_nonref in level]
        plot_levels_nonref_2 = [level for level in plot_levels if substring_nonref_2 in level]
        plot_levels_ordered = plot_levels_ref_2 + plot_levels_ref_nonref + plot_levels_nonref_2
        fig, ax = plt.subplots()
        violin_df[gene] = violin_df[gene].sparse.to_dense()
        sns.violinplot(x=violin_df[marker], y=violin_df[gene], order=order, ax=ax)
        sns.stripplot(x=violin_df[marker], y=violin_df[gene], color='black', size=5, order=order, ax=ax, jitter=True)
        ax.set_xticklabels(plot_levels_ordered)
        ## Add text
        beta=row['BETA']
        fdr=row['bonf_qvalue_across_genes']
        text = f"beta = {float(beta):.2f}, FDR={fdr:.2e}"
        ax.text(1, 1.1 * ax.get_ylim()[1], text, ha='center', va='bottom')
        # Add options 
        ax.set_xlabel(marker)
        ax.set_ylabel(gene)
        # save
        plt.savefig(f"{outdir}/{marker}_{gene}_trans_sns.png")
        plt.clf()
    except Exception as e:
        # Print or log the error message
        print(f"Error processing row {r}: {e}")


# Have a look at the number of trans-eGenes per variant
var_counts = res_sig['MarkerID'].value_counts()

# Plot a histogram of this
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
plt.hist(var_counts, bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
# Add labels and title
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Trans-eGenes per cis-eQTL')
plt.savefig(f"{outdir}/sig_trans_eGenes_per_variant.png")
plt.clf()

# Plot the number of eQTLs/gene and the relationship with gene expression - is it driven by lowly expressed genes? Correlate with MAF too
expression_df = pd.DataFrame.sparse.from_spmatrix(adata.X)
expression_df.columns=adata.var.index.values
expression_df.index = adata.obs.index.values
subexpr = expression_df.iloc[:,expression_df.columns.isin(res.Gene.values)]
subexprmeans = subexpr.mean()
subexprmeans = pd.DataFrame({'Gene': subexprmeans.index, 'Mean': subexprmeans.values})
res = res.merge(subexprmeans, on="Gene")
res['fdr_less_0.01'] = np.where(res['bonf_qvalue_across_genes'] < 0.01, 'Yes', 'No')
res['MAF'] = np.where(res['AF_Allele2'] > 0.5, 1 - res['AF_Allele2'], res['AF_Allele2'])
res['logMean'] = np.log10(res['Mean'])

plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='MAF', y='logMean', data=res, hue="fdr_less_0.01")
plt.xlabel('MAF')
plt.ylabel('log10(Mean expression of gene)')
plt.xlim(0, 0.5)
plt.title(f"Detection of trans-eQTLs (MAF vs gene expression): FDR<0.01")
plt.savefig(f"{outdir}/trans_detection_MAF_expr_py.png")
plt.clf()

# Plot a distribution of MAF across significance thresholds
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
groups = np.unique(res['fdr_less_0.01'])
for g in groups:
    sns.distplot(res[res['fdr_less_0.01'] == g].MAF, hist=False, rug=True, label=g, kde_kws={'linewidth': 1.5})

plt.legend(bbox_to_anchor=(1.0, 1.0))
plt.title('Distribution of trans-eQTL MAF across significance (FDR<0.01)')
plt.xlabel('MAF')
plt.xlim(0,0.5)
#plt.axvline(x = 0.7, color = 'red', linestyle = '--', alpha = 0.5)
plt.savefig(f"{outdir}/dist_trans_MAF_significance_py.png", bbox_inches='tight')
plt.clf()

# How does beta value vary with the MAF and gene expression? 
res['absBETA'] = np.abs(res['BETA'])
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='MAF', y='logMean', data=res, hue="absBETA", palette="viridis", alpha=0.8)
plt.xlabel('MAF')
plt.ylabel('log10(Mean expression of gene)')
plt.xlim(0, 0.5)
plt.title(f"Detection of trans-eQTLs (MAF vs gene expression), colored by BETA")
plt.savefig(f"{outdir}/trans_detection_MAF_expr_beta_py.png")
plt.clf()

# BETA vs significance
res['log10_p.value'] = -np.log10(res['p.value'].astype('float'))
res['log10_bonf_qvalue_across_genes'] = -np.log10(res['bonf_qvalue_across_genes'])
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='Mean', y='log10_p.value', data=res)
plt.xlabel('logMean')
plt.ylabel('log10(p-value)')
plt.xlim(0, 0.5)
plt.savefig(f"{outdir}/trans_detection_significance_expr_py.png")
plt.clf()

# Plot the detection of trans vs sparsity (number of expressing cells and samples)
counts = expression_df
counts = counts.loc[:, counts.columns.isin(res['Gene'])]
counts.astype(bool).sum(axis=0)
counts_per_sample = counts.astype(bool).sum(axis=0)
counts_per_sample = pd.DataFrame(counts_per_sample)
counts_per_sample = counts_per_sample.reset_index()
counts_per_sample.columns=["Gene", "n_expr_cells"]
res = res.merge(counts_per_sample, on="Gene")
res['gene_expr_sparsity_cell'] = res['n_expr_cells'].astype('float') / counts.shape[0]
res['log10_p.value'] = -np.log10(res['p.value'].astype('float'))
plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='gene_expr_sparsity_cell', y='log10_p.value', data=res)
plt.xlabel('Proportion of expressing cells')
plt.ylabel('log10(p-value)')
plt.xlim(0,1)
plt.savefig(f"{outdir}/trans_detection_significance_expr_sparsity_cell_py.png")
plt.clf()

# Find the absolute number of homozygous non ref (2) samples for each genotype. Plot this against the significance/beta
G = read_plink1_bin(f"{geno_dir}/plink_genotypes_cis_{level}.bed")
G_id = pd.read_csv(f"{geno_dir}/plink_genotypes_cis_{level}.bim", sep = "\t", header=None)
G['variant'] = G_id.iloc[:,1]
G = G[G.sample.isin(adata.obs.Corrected_genotyping_ID)]
G = G[:,G['variant'].isin(res.MarkerID)]
# Plot the expression of these genes
genotype_df = G.to_dataframe(name='genotype').reset_index()
genotype_df = genotype_df.pivot(index='sample', columns='variant', values='genotype')
expression_df = pd.DataFrame.sparse.from_spmatrix(adata.X)
expression_df.columns=adata.var.index.values
expression_df.index = adata.obs.index.values
expression_df = expression_df.iloc[:,expression_df.columns.isin(res_plot.Gene)]
genotype_df = genotype_df.reset_index()
genotype_df = genotype_df.rename(columns={'sample': 'Corrected_genotyping_ID'})
nonref_dict = {}
for column in genotype_df.columns:
    if column.startswith('chr'):
        nonref_dict[column] = genotype_df[column].astype(str).eq('0.0').sum()

count_hom_nonref = pd.DataFrame(list(nonref_dict.items()), columns=['MarkerID', 'Count_hom_nonref'])
res = res.merge(count_hom_nonref, on="MarkerID")

plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='Count_hom_nonref', y='log10_p.value', data=res)
plt.xlabel('Number of homozygous non-ref samples')
plt.ylabel('log10(p-value)')
plt.savefig(f"{outdir}/trans_detection_significance_hom_non_ref_py.png")
plt.clf()

plt.figure(figsize=(8, 6))
fig,ax = plt.subplots(figsize=(8,6))
sns.scatterplot(x='Count_hom_nonref', y='absBETA', data=res)
plt.xlabel('Number of homozygous non-ref samples')
plt.ylabel('absolute(BETA)')
plt.savefig(f"{outdir}/trans_detection_beta_hom_non_ref_py.png")
plt.clf()
