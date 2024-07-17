#### Bradley Jan 2024
#### Summary of the SAIGE vs Tensor eQTL effects across multiple cell-types

# Load packages
# Import libraries
import os
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
import argparse
from sklearn.preprocessing import LabelEncoder
from scipy.stats import entropy
from plotnine import ggplot, aes, geom_point, labs, theme, element_blank, theme_minimal
print("Loaded libraries")


# Define the resolutions and annotations to work with
chromosomes = "1-22" # Will be something like "1-5", so divide this by '-'
start, end = map(int, chromosomes.split('-'))
chr = range(start, end + 1)
formatted_range = chromosomes
aggregate_on = "label__machine"
levels=["B_cell_activated", "Dendritic_cell", "Mac_resident_IL10RAplus",  "T_cell_CD4_CD40LGplus_2", "T_cell_CD8_1", "Tuft_cell", "Stem_cell_LGR5plus", "Enterocyte_precursor_crypt_OLFM4plus_KRT20plusplus"]
basedir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"

# Load anndata to get the number of cells tested per level
adata = ad.read_h5ad(phenotype__file)
adata = adata[adata.obs[aggregate_on].isin(levels)]
level_counts = pd.DataFrame(np.unique(adata.obs[aggregate_on], return_counts=True)).T
level_counts.columns = ["Level", "nCells"]
label_cat = adata.obs[['label__machine', 'category__machine']]
label_cat = label_cat.reset_index()
label_cat = label_cat[['label__machine', 'category__machine']]
label_cat.drop_duplicates(inplace=True)
label_cat = label_cat.rename(columns={aggregate_on: 'Level'})

# For each annotation, load the results for the desired chromosomes
print("Loading in the sc vs pb results for the chromosome choice")
multi_ct=[]
for l in levels:
    print(f"Reading in data for {l}")
    saige_tensor_temp = pd.read_csv(f"{basedir}/{aggregate_on}/{l}/summary_tables/min_qvalue_per_common_gene_saige_tensor_{formatted_range}.txt", sep = "\t")
    saige_tensor_temp['label'] = l
    columns_to_convert = ['qvalue_across_genes_saige', 'qvalue_across_genes_tensor', 'p.value', 'pval_nominal']
    saige_tensor_temp[columns_to_convert] = saige_tensor_temp[columns_to_convert].astype(float)
    multi_ct.append(saige_tensor_temp)

# Count the significance (FDR across genes < 0.05)
fdr_df = []
# Iterate over the list of DataFrames
for i, df in enumerate(multi_ct):
    # Count the number of rows where the conditions are met
    saige_count = (df['qvalue_across_genes_saige'] < 0.05).sum()
    tensor_count = (df['qvalue_across_genes_tensor'] < 0.05).sum()
    fc = saige_count / tensor_count if tensor_count != 0 else float('nan')
    # Append the counts to the DataFrame
    fdr_df.append({'Level': levels[i], 'Saige_Count': saige_count, 'Tensor_Count': tensor_count, 'nSaige_nTensor': fc})

fdr_df = pd.DataFrame(fdr_df)
fdr_df = fdr_df.merge(level_counts, on="Level")
fdr_df = fdr_df.merge(label_cat, on="Level")
print(fdr_df)

# Calculate the shannon entropy of each cluster and add this to the results df
dense_matrix = adata.X.toarray()
adata.obs['shannon_entropy'] = np.apply_along_axis(lambda x: entropy(x, base=2), axis=1, arr=dense_matrix)
# Calculate min/mean per sample
grouped_stats = adata.obs.groupby([aggregate_on, 'sanger_sample_id'])['shannon_entropy'].agg(['mean', 'median', 'min', 'max', lambda x: np.ptp(x)]).reset_index()
grouped_stats['range'] = grouped_stats['max'] - grouped_stats['min']
result_df = grouped_stats.pivot(index='label__machine', columns='sanger_sample_id', values=['mean', 'median', 'min', 'max', 'range']).reset_index()
result_df = grouped_stats.groupby('label__machine').agg({
    'mean': 'mean',
    'median': 'mean',
    'min': 'mean',
    'max': 'mean',
    'range': 'mean'
}).reset_index()

# Rename the columns for clarity
result_df.columns = result_df.columns[0:1].tolist() + ['shannon_' + col for col in result_df.columns[1:]]
result_df = result_df.rename(columns={aggregate_on: 'Level'})
fdr_df = fdr_df.merge(result_df, on="Level")

# Also caculate the average distances between cells/cell-type
from sklearn.metrics import pairwise_distances
means = []
for l in levels:
    print(l)
    subset_expression = adata.X[adata.obs[aggregate_on]==l]
    distances = pairwise_distances(subset_expression, metric='euclidean')
    average_distance = np.mean(distances)
    means.append(average_distance)

euclid = pd.DataFrame({'Level': levels, 'average_euclid_distance': means})
fdr_df = fdr_df.merge(euclid, on="Level")

# Plot the relationship between number of genes detected by each method and the number of cells/shannon entropy
outdir=f"{basedir}/overall_plots"
if os.path.exists(outdir) == False:
    os.makedirs(outdir, exist_ok=True)


plt.figure(figsize=(10, 6))
sns.scatterplot(data=fdr_df, x='nCells', y='Saige_Count', hue='category__machine', marker='o', s=150, edgecolor='black')
sns.scatterplot(data=fdr_df, x='nCells', y='Tensor_Count', marker='o', s=150, edgecolor='black', facecolor='none')
# Set labels and title
plt.xlabel('nCells')
plt.ylabel('n eGenes (FDR < 0.05)')
# Annotate the maximum point for each level
for idx, row in fdr_df.loc[fdr_df.groupby('Level')['Saige_Count'].idxmax()].iterrows():
    level = row['Level']
    max_fc = row['nSaige_nTensor']
    max_x = row['nCells']
    max_y = max(row['Saige_Count'], row['Tensor_Count'])
    print(f'{max_fc:.2f}')
    plt.text(max_x, max_y, f'Level:{level}\nFC:{max_fc:.2f}', ha='left', va='bottom', fontsize=8, color='black')

plt.title('SAIGEQTL vs TensorQTL')
plt.savefig(outdir + '/n_eGenes_SAIGE_Tensor_ncells.png', bbox_inches='tight')
plt.clf()

# Relationship with shannon entropy
sns.scatterplot(data=fdr_df, x='shannon_range', y='Saige_Count', hue='category__machine', marker='o', s=150, edgecolor='black')
sns.scatterplot(data=fdr_df, x='shannon_range', y='Tensor_Count', marker='o', s=150, edgecolor='black', facecolor='none')
# Set labels and title
plt.xlabel('Mean Shannon Entropy')
plt.ylabel('n eGenes (FDR < 0.05)')
plt.title('SAIGEQTL vs TensorQTL')
plt.savefig(outdir + '/n_eGenes_SAIGE_Tensor_shannon.png', bbox_inches='tight')
plt.clf()

# Relationship with shannon entropy
sns.scatterplot(data=fdr_df, x='average_euclid_distance', y='Saige_Count', hue='category__machine', marker='o', s=150, edgecolor='black')
sns.scatterplot(data=fdr_df, x='average_euclid_distance', y='Tensor_Count', marker='o', s=150, edgecolor='black', facecolor='none')
# Set labels and title
plt.xlabel('Average cell-cell euclidean distance')
plt.ylabel('n eGenes (FDR < 0.05)')
plt.title('SAIGEQTL vs TensorQTL')
plt.savefig(outdir + '/n_eGenes_SAIGE_Tensor_euclid.png', bbox_inches='tight')
plt.clf()
