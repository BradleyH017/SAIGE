# Bradley Jan 2024
# Checking the conditional analysis results wrt p-value correction

# Load libraries
import matplotlib as mp
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import scanpy as sc
import statsmodels.stats.multitest as smt
from scipy.stats import linregress
from scipy import stats
import scipy.stats as stats
import statsmodels.api as sm 
import seaborn as sns
import matplotlib.pyplot as plt

# Define options for checking
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
aggregate_on="label__machine"
level="T_cell_CD8_1"
outdir=f"{general_file_dir}/{aggregate_on}/{level}/summary_plots"
n_expr_pcs=15

methods = ["bh"]
sigres = []
for m in methods:
    # Load the conditional results
    sig = pd.read_csv(f"{general_file_dir}/{aggregate_on}/{level}/conditional/{m}/all_conditionally_independent_effects_q_less_0.05.txt", sep = "\t", header=None)
    # Get header from other file
    temp = pd.read_csv(f"{general_file_dir}/{aggregate_on}/{level}/conditional/{m}/all_conditionally_independent_effects.txt", sep = "\t")
    sig.columns = temp.columns
    sigres.append(sig)
    # Plot a histogram of effects / n effects per gene
    count_sig_per_gene = pd.DataFrame(np.unique(sig['Gene'], return_counts=True)).T
    count_sig_per_gene.columns = ["Gene", "nSignals"]
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.hist(count_sig_per_gene['nSignals'], bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    # Add labels and title
    plt.xlabel('nSignals')
    plt.ylabel('Frequency')
    plt.title('nSignals per gene')
    plt.savefig(f"{outdir}/nSignals_per_gene_{m}.png")
    plt.clf()
    # Plot violings of the nominal p-values per variants
    sig['log10_p_nominal'] = -np.log10(sig['p.value'])
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='round', y='log10_p_nominal', data=sig, split=True, inner='quartile')
    # Add lines connecting points within the same 'Gene'
    sns.pointplot(x='round', y='log10_p_nominal', data=sig, dodge=0.2, markers='o', color='grey', markersize=4, alpha=0.7, hue='Gene', legend=False)
    plt.title('Nominal p-value per independent effect')
    plt.savefig(f"{outdir}/nom_p_per_signal_per_gene_{m}.png")
    plt.clf()

# Plot the LD between these variants
# Calculate LD:
# conda activate bcf
# cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/genotypes
# plink --bfile plink_genotypes_cis_T_cell_CD8_1 --r2 --out ld_cis_T_cell_CD8_1 --tab --ld-window-r2 0 --ld-window-kb 100000000000 --ld-window 1000000

# Load this in
ld = pd.read_csv(f"{general_file_dir}/genotypes/ld_cis_{level}.ld", delim_whitespace=True)

res_ld=[]
for index, m in enumerate(methods):
    # Add this to the significant results
    sig_secondaries = sigres[index][sigres[index]['round'].astype('float') > 1]
    # Subset for those with LD results. NOTE: This should not be neccessary as LD should be calculated for all variants
    sig_secondaries = sig_secondaries[sig_secondaries['MarkerID'].isin(ld['SNP_A']) | sig_secondaries['MarkerID'].isin(ld['SNP_B'])]
    signal_lds = []
    genes = np.unique(sig_secondaries['Gene'])
    for g in genes:
        print(g)
        res = sig_secondaries[sig_secondaries['Gene'] == g]
        dummy = pd.DataFrame(np.empty((res.shape[0], max(res['round']-1)), dtype=object))
        dummy.index = [f'{i}' for i in range(1, res.shape[0]+1)]
        dummy.columns = [f'{i}' for i in range(2, res.shape[0]+2)]
        for ref_round in dummy.index:
            for query_round in dummy.columns:
                varref = sig[(sig['round'] == int(ref_round)) & (sig['Gene'] == g) ]['MarkerID'].values[0]
                varquery = sig[(sig['round'] == int(query_round)) & (sig['Gene'] == g)]['MarkerID'].values[0]
                if varref != varquery:
                    ld_to_add = ld[(ld['SNP_A'].values == varquery) & (ld['SNP_B'].values == varref)]['R2'].values
                    # If this doesn't find it
                    if len(ld_to_add) == 0:
                        ld_to_add = ld[(ld['SNP_B'].values == varquery) & (ld['SNP_A'].values == varref)]['R2'].values
                    # Add to df
                    dummy.iloc[int(ref_round)-1, int(query_round)-2] = ld_to_add
        dummy['ref_round'] = dummy.index
        dummy = pd.melt(dummy, id_vars='ref_round', value_name="r2")
        dummy = dummy.rename(columns={'variable': 'query_round'})
        dummy['Gene'] = g
        signal_lds.append(dummy)
    all_lds = pd.concat(signal_lds)
    all_lds = all_lds[all_lds['r2'].notna()]
    res_ld.append(all_lds)
    # Now extract the mean, maximum values for each comparison
    mean_combos = pd.DataFrame(np.empty((4, 4), dtype=object))
    mean_combos.index = [f'{i}' for i in range(1, 5)]
    mean_combos.columns = [f'{i}' for i in range(2, 6)]
    max_combos = pd.DataFrame(np.empty((4, 4), dtype=object))
    max_combos.index = [f'{i}' for i in range(1, 5)]
    max_combos.columns = [f'{i}' for i in range(2, 6)]
    for r in mean_combos.index:
        for c in mean_combos.columns:
            if int(c) > int(r):
                print(f"round{r} vs round{c}")
                temp = all_lds[(all_lds['ref_round'] == r) & (all_lds['query_round'] == c)]['r2']
                mean_combos.loc[r,c] = sum(temp) / len(temp)
                max_combos.loc[r,c] = max(temp)
    # Plot each of these
    plt.figure(figsize=(8, 8))
    mean_combos = mean_combos.apply(pd.to_numeric, errors='coerce')
    sns.heatmap(mean_combos, annot=True, cmap='viridis', fmt='.2f', square=True, vmin=0, vmax=1)
    # Set axis labels and title
    plt.xlabel('Round')
    plt.ylabel('Round')
    plt.title('Mean r2')
    plt.savefig(f"{outdir}/mean_ld_independent_pairs_{m}.png")
    plt.clf()
    # Also plot max
    max_combos = pd.DataFrame(max_combos)
    plt.figure(figsize=(8, 8))
    max_combos_numeric = max_combos.applymap(lambda x: float(x) if x is not None else x)
    sns.heatmap(max_combos_numeric, annot=True, cmap='viridis', fmt='.2f', square=True, vmin=0, vmax=1)
    # Set axis labels and title
    plt.xlabel('Round')
    plt.ylabel('Round')
    plt.title('Max r2')
    plt.savefig(f"{outdir}/max_ld_independent_pairs_{m}.png")
    plt.clf()

# Make a manhattan of some good and bad loci per method
# Sum r2 per gene - pick 3 high LD and 3 low LD
test_genes = np.unique(res_ld[0][res_ld[0]['query_round'].astype(float) > 3]['Gene'])
sum_df = res_ld[0][res_ld[0]['Gene'].isin(test_genes)].groupby('Gene')['r2'].sum().reset_index()
sum_df = sum_df.sort_values(by='r2', ascending=False)
high = sum_df.head(3)['Gene'].values
low = sum_df.tail(3)['Gene'].values
plot_genes = set(high).union(set(low))
# Plot a manhattan for each of these, highlighting the conditionally independent variants
for g in plot_genes:
    print(f"Plotting {g}")
    # Load the results for this chromosome
    gene_coords_file = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
    gene_coords = pd.read_csv(gene_coords_file, sep = "\t")
    columns_to_read = list(range(19))
    chr = np.unique(sigres[0][sigres[0]['Gene'] == g]['CHR'])[0]
    res_chr = pd.read_csv(f"{general_file_dir}/{aggregate_on}/{level}/cis/chr{str(chr)}_nPC_{str(n_expr_pcs)}.txt.gz", compression='gzip', delimiter='\t', usecols=columns_to_read)
    gene_res = res_chr[res_chr['Gene'] == g]
    sigres_gene = sigres[0][sigres[0]['Gene'] == g]
    marker_round_map = dict(zip(sigres_gene['MarkerID'], sigres_gene['round']))
    gene_res['round'] = gene_res['MarkerID'].map(marker_round_map).fillna('No')
    # Plot a manhattan of these
    gene_res['color'] = gene_res['round'].apply(lambda x: 'navy' if x == 'No' else 'red')
    gene_res['log10p'] = -np.log10(gene_res['p.value'].astype(float))
    gene_coords_window=gene_coords[(gene_coords['start'] > min(gene_res['POS'])) & (gene_coords['end'] < max(gene_res['POS'])) & (gene_coords['chromosome'] == str(chr)) & gene_coords['feature_id'].isin(res_chr['Gene'])]
    #
    plt.figure(figsize=(10, 6))
    fig, ax1 = plt.subplots(figsize=(10, 6))
    sns.scatterplot(data=gene_res, x='POS', y='log10p', hue='color', palette={'navy': 'navy', 'red': 'red'}, s=50)
    # Annotate points with 'round' values
    for index, row in gene_res[gene_res['round'] != 'No'].iterrows():
        plt.annotate(row['round'], (row['POS'], row['log10p']), color='black', fontsize=10, ha='center', va='bottom')
    #
    plt.xlabel('position')
    plt.ylabel('-log10(p-value)')
    plt.title(f"Conditionally independent variants - {g}")
    plt.ylim([0,max(gene_res['log10p'])+2])
    #
    if g in high:
        group="low"
    else:
        group="high"
    #
    plt.savefig(f"{outdir}/{g}_conditional_variants_group_{group}_LD.png")
    plt.clf()
    #
    
    
    
# Also plot the lds
ld_test = all_lds[all_lds['Gene'] == g]
ld_test = ld_test.pivot(index='ref_round', columns='query_round', values='r2')
plt.figure(figsize=(8, 8))
ld_test = ld_test.apply(pd.to_numeric, errors='coerce')
sns.heatmap(ld_test, annot=True, cmap='viridis', fmt='.2f', square=True, vmin=0, vmax=1)
# Set axis labels and title
plt.xlabel('Round')
plt.ylabel('Round')
plt.title('Mean r2')
plt.savefig(f"{outdir}/{g}_ld_independent_pairs.png")
plt.clf()


# Plot the relationship between the number of signals detected and the minimum p value
