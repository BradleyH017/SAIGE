#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-23'
__version__ = '0.0.1'

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
print("Loaded libraries")

# Define genotype replacement functions
def replace_with_ref_nonref(value, ref, nonref):
    if value == '0.0':
        return nonref * 2  # Repeat 'nonref' twice for '0'
    elif value == '1.0':
        return ref + nonref  # Paste 'ref' and 'nonref' for '1'
    elif value == '2.0':
        return ref * 2  # Repeat 'ref' twice for '2'
    return value

# Inherit option
def parse_options():    
    # Inherit options
    parser = argparse.ArgumentParser(
            description="""
                Prepping files for the SAIGEQTL analysis
                """
        )

    parser.add_argument(
            '-c', '--catdir',
            action='store',
            dest='catdir',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-ch', '--chromosomes_cis',
            action='store',
            dest='chromosomes',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-pbd', '--pb_dir',
            action='store',
            dest='pb_dir',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-ph', '--phenotype__file',
            action='store',
            dest='phenotype__file',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-o', '--general_file_dir',
            action='store',
            dest='general_file_dir',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-nx', '--n_expr_pcs',
            action='store',
            dest='n_expr_pcs',
            required=True,
            help=''
        )
    
    return parser.parse_args()

def main():
    print("Inheriting running options")
    inherited_options = parse_options()
    catdir = inherited_options.catdir
    # catdir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles/label__machine/T_cell_CD8_1"
    chromosomes = inherited_options.chromosomes # Will be something like "1-5", so divide this by '-'
    start, end = map(int, chromosomes.split('-'))
    chr = range(start, end + 1)
    formatted_range = chromosomes
    # Extract the other options from the catdir
    # Specify sc/pseudo-bulk files for this chromosome
    options = catdir.split("/")
    aggregate_on= [element for element in options if "__machine" in element][0]
    level = options[-1]
    general_file_dir=inherited_options.general_file_dir
    #general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
    geno_dir=f"{general_file_dir}/genotypes"
    pb_dir=inherited_options.pb_dir
    # pb_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/Pseudobulk-TensorQTL/2023_09_28-TI_fr003-plate123"
    n_expr_pcs=inherited_options.n_expr_pcs
    # n_expr_pcs="5"
    phenotype_file=inherited_options.phenotype__file
    # phenotype_file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"

    # Print the options for running this script
    import logging
    logging.info(f"catdir: {catdir}")
    logging.info(f"chromosomes: {chromosomes}")
    logging.info(f"general_file_dir: {general_file_dir}")
    logging.info(f"pb_dir: {pb_dir}")
    logging.info(f"phenotype__file: {phenotype_file}")
    logging.info(f"n_expr_pcs: {n_expr_pcs}")

    # Prep the outdir if not already
    outdir=f"{catdir}/summary_plots"
    if os.path.exists(outdir) == False:
        os.makedirs(outdir, exist_ok=True)

    tabdir=f"{catdir}/summary_tables"
    if os.path.exists(tabdir) == False:
        os.makedirs(tabdir, exist_ok=True)

    # Start with a summary of the cis-eQTL analysis for the choromosomes specified
    print("~~~~~~~~~~~~~~~~~~ Summarising the cis-eQTL analysis ~~~~~~~~~~~~~~~~~~")
    print("Loading in the results for the chromosome choice")
    # Load in the sc_res
    sc_res = []
    for c in chr:
        print(c)
        fname = f"{catdir}/cis/chr{str(c)}_nPC_{n_expr_pcs}.txt.gz"
        f = pd.read_csv(fname, sep = "\t", usecols=range(17), compression='gzip')
        sc_res.append(f)

    sc_res = pd.concat(sc_res, axis=0, ignore_index=True)
    sc_res = sc_res.dropna(subset=['Gene'])
    mask = sc_res['Gene'].notna()
    sc_res = sc_res[mask & sc_res['Gene'].str.contains("ENSG")]
    sc_res['MarkerID'] = sc_res['MarkerID'].astype(str).replace('nan', '')
    sc_res['Gene'] = sc_res['Gene'].astype(str).replace('nan', '')

    # Combine the phenotype/variant columns
    sc_res['variant_phenotype'] = sc_res['MarkerID'] + "_" + sc_res['Gene']
    sc_res = sc_res.drop_duplicates(subset=['variant_phenotype'])

    # Plot quality control plots for the single cell expression only
    print("Plotting p-value histogram")
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

    print("Plotting QQ")
    # Also QQ plot
    plt.figure(figsize=(8, 6))
    obs = sc_res['p.value']
    n = len(obs)
    obs_p = -np.log10(np.sort(obs))
    th_p = np.arange(1/float(n), 1 + (1/float(n)), 1/float(n))
    th_p = -np.log10(th_p)
    if len(th_p) > len(obs_p):
        # Truncate th_p to the same length as obs_p
        th_p = th_p[:len(obs_p)]
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

    print("Plotting diagnostic plots of results")
    print("Relationship between absolute effect size and significance per gene")
    # Plot the relationship between absolute beta and variant MAF (top hit per gene)
    sc_res['log10_pval'] = -np.log10(sc_res['p.value'].astype('float'))
    sc_res['abs_beta'] = np.abs(sc_res['BETA'].astype('float'))
    sc_res['AF_Allele2'] = sc_res['AF_Allele2'].astype('float')
    sc_res['MAF'] = np.where(sc_res['AF_Allele2'] > 0.5, 1 - sc_res['AF_Allele2'], sc_res['AF_Allele2'])
    max_indices = sc_res.groupby('Gene')['log10_pval'].idxmax()
    sc_res_top_per_gene = sc_res.loc[max_indices]
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='MAF', y='abs_beta', data=sc_res_top_per_gene)
    plt.xlabel('MAF')
    plt.ylabel('Absolute effect size')
    plt.xlim(0, 0.5)
    plt.title(f"Chromosome {formatted_range} - top hit/gene")
    plt.savefig(f"{outdir}/absolute_effect_vs_maf_saige_chr{formatted_range}_py.png")
    plt.clf()

    # Number of variants per gene
    print("The number of variants being tested per gene")
    var_count = sc_res['Gene'].value_counts()
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.hist(var_count, bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    # Add labels and title
    plt.xlabel('Number of variants tested per gene')
    plt.ylabel('Frequency')
    plt.title('Histogram of variants tested per gene')
    plt.savefig(f"{outdir}/variants_per_gene{formatted_range}_py.png")
    plt.clf()

    # Number of genes per variant (non0)
    print("Number of genes being tested per variant (and the number of siginificant hits per variant)")
    gene_count = sc_res['MarkerID'].value_counts()
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.hist(gene_count, bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    # Add labels and title
    plt.xlabel('Number of genes tested per variant')
    plt.ylabel('Frequency')
    plt.title('Histogram of genes tested per variant')
    plt.savefig(f"{outdir}/genes_per_variant_{formatted_range}_py.png")
    plt.clf()

    # Number of significant hits per variant (vs number of genes tested)
    print("The number of significant (within gene q<0.05) hits per variant with effect")
    # Do a quick q-value correction within genes
    def apply_multiple_testing_correction(df, group_column, pvalue_column):
        # Group by the specified column
        grouped_df = df.groupby(group_column)
        # Define a function to apply multiple testing correction for each group
        def correction(group):
            _, q_values, _, _ = smt.multipletests(group[pvalue_column], method='fdr_bh')
            group['q_value_within_gene'] = q_values
            return group
        # Apply the correction function to each group
        df = grouped_df.apply(correction)    
        return df

    sc_res = apply_multiple_testing_correction(sc_res, 'Gene', 'p.value')
    sig_var_count = sc_res[sc_res['q_value_within_gene'] < 0.05]['MarkerID'].value_counts()
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.hist(sig_var_count, bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    # Add labels and title
    plt.xlabel('Number of hits per variant')
    plt.ylabel('Frequency')
    plt.title('Histogram of significant effects per variant (within gene q<0.05)')
    plt.savefig(f"{outdir}/sig_genes_per_sig_variant_{formatted_range}_py.png")
    plt.clf()

    ## Number of signals per gene
    #print("The number of independent signals per gene")
    #conditional = pd.read_csv(f"{catdir}/conditional/all_conditionally_independent_effects_q_less_0.05.txt", sep = "\t", index_col=None, header=None)
    #condtional=conditional.iloc[:,:-2]
    #names=['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'AC_Allele2', 'AF_Allele2', 'MissingRate', 'BETA', 'SE', 'Tstat', 'var', 'p.value', 'p.value.NA', 'Is.SPA', 'BETA_c', 'SE_c', 'Tstat_c', 'var_c', 'p.value_c', 'p.value.NA_c', 'N', 'qvalues_c', 'lfdr_c', 'round', 'Gene']
    #conditional.set_index('CHR', inplace=True)
    #conditional = conditional[['POS'] + list(conditional.columns[:-1])]

    #conditional = conditional[conditional['qvalues_c'] < 0.05]
    #conditional_gene_counts = conditional['Gene'].value_counts()
    #plt.figure(figsize=(8, 6))
    #fig,ax = plt.subplots(figsize=(8,6))
    #plt.hist(conditional_gene_counts, bins=40, color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    ## Add labels and title
    #plt.xlabel('Number of independent signals hits per gene')
    #plt.ylabel('Frequency')
    #plt.title('Histogram of conditional hits per gene')
    #plt.savefig(f"{outdir}/independent_effects_per_gene_{formatted_range}_py.png")
    #plt.clf()

    print("~~~~~~~~~~~~~~~~~~ Comparison of cis-eQTLs between SAIGE and TENSOR ~~~~~~~~~~~~~~~~~~")
    # Now compare this with the pseudo-bulk results
    print("Loading the pseudobulk results for these chromosomes")
    pb_res = []
    pb_base=f"{pb_dir}/{aggregate_on}-{level}-dMean"
    for c in chr:
        print(c)
        fname = f"{pb_base}/cis_nominal1.cis_qtl_pairs.chr{str(c)}.tsv"
        f = pd.read_csv(fname, sep = "\t")
        pb_res.append(f)

    pb_res = pd.concat(pb_res, axis=0, ignore_index=True)
    pb_res['variant_phenotype'] = pb_res['variant_id'] + "_" + pb_res['phenotype_id']
        
    # Before filtering, FDR correct these reslts within and across genes
    pb_res['qvalue_within_gene'] = pb_res.groupby('phenotype_id')['pval_nominal'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
    pb_res_min_qvalue_rows = pb_res.loc[pb_res.groupby('phenotype_id')['qvalue_within_gene'].idxmin()]
    pb_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(pb_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]
    sc_res['p.value'] = sc_res['p.value'].astype('float')
    sc_res['qvalue_within_gene'] = sc_res.groupby('Gene')['p.value'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
    sc_res_min_qvalue_rows = sc_res.loc[sc_res.groupby('Gene')['qvalue_within_gene'].idxmin()]
    sc_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(sc_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]
    nsaige_sig = (sc_res_min_qvalue_rows['sc_res_min_qvalue_rows'] < 0.05).sum()
    ntensor_sig = (sc_res_min_qvalue_rows['qvalue_across_genes'] < 0.05).sum()

    # Plot the overlap of these results
    saige_hits = set(sc_res_min_qvalue_rows[sc_res_min_qvalue_rows['qvalue_across_genes'] < 0.05]['Gene'].tolist())
    tensor_hits = set(pb_res_min_qvalue_rows[pb_res_min_qvalue_rows['qvalue_across_genes'] < 0.05]['phenotype_id'].tolist())
    venn_labels = {'100': 'SAIGE-QTL', '010': 'TensorQTL', '110': 'Intersection'}
    plt.figure(figsize=(8, 6))
    venn2(subsets=[saige_hits, tensor_hits], set_labels=('SAIGE-QTL', 'TensorQTL'), set_colors=('skyblue', 'lightgreen'), alpha=0.7)
    plt.savefig(f"{outdir}/venn_saige_tensor_fdr_0.05_across_genes_all_not_intersection_chr{formatted_range}_py.png")
    plt.clf()

    # Save these datasets
    sc_res_min_qvalue_rows.to_csv(f"{tabdir}/min_qvalue_per_gene_saige_all_tests_{formatted_range}.txt", sep = "\t", index=False) 
    pb_res_min_qvalue_rows.to_csv(f"{tabdir}/min_qvalue_per_gene_saige_all_tests_{formatted_range}.txt", sep = "\t", index=False) 


    ########### Now subset to the common tests to directly compare the models ###########
    # Find common tests
    common = set(pb_res['variant_phenotype']).intersection(sc_res['variant_phenotype'])

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
    sns.scatterplot(x='slope', y='BETA', data=merged)
    plt.plot([-max_effect, max_effect], [-max_effect, max_effect], color='red', linestyle='dashed', label='x=y Line')
    plt.xlabel('Effect size - TensorQTL')
    plt.ylabel('Effect size - SAIGE')
    plt.xlim(-max_effect, max_effect)
    plt.ylim(-max_effect, max_effect)
    plt.text(max_effect, -max_effect, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
    plt.title(f"Chromosome {formatted_range} only")
    # Save the plot
    plt.savefig(f"{outdir}/cor_betas_pb_chr{formatted_range}_py.png")
    plt.clf()

    # Do the same but now look at p-value
    max_effect = max([np.abs(max(-np.log10(merged['pval_nominal']))), np.abs(max(-np.log10(merged['p.value'].astype('float'))))])
    merged['tensor_nlog10'] = -np.log10(merged['pval_nominal'])
    merged['saige_nlog10'] = -np.log10(merged['p.value'].astype(float))
    rho, p_value = pearsonr(merged['saige_nlog10'], merged['tensor_nlog10'])
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='tensor_nlog10', y='saige_nlog10', data=merged)
    plt.plot([0, max_effect], [0, max_effect], color='red', linestyle='dashed', label='x=y Line')
    plt.xlabel('-log10(TensorQTL p-value)')
    plt.ylabel('-log10(SAIGEQTL p-value)')
    plt.xlim(0, max_effect)
    plt.ylim(0, max_effect)
    plt.title(f"Chromosome {formatted_range} only")
    plt.text(max_effect, 0, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
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


    # Is the win rate dependent on MAF
    # Compare the absolute difference in beta across MAF
    merged['AF_Allele2'] = merged['AF_Allele2'].astype('float')
    merged['absolute_diff_betas'] = np.abs(merged['slope'] - merged['BETA'])
    merged['MAF'] = np.where(merged['AF_Allele2'] > 0.5, 1 - merged['AF_Allele2'], merged['AF_Allele2'])
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='MAF', y='absolute_diff_betas', data=merged)
    plt.xlabel('MAF')
    plt.ylabel('Absolute difference in effect size (SAIGE-QTL vs TensorQTL)')
    plt.xlim(0, 0.5)
    plt.title(f"Chromosome {formatted_range} only")
    plt.savefig(f"{outdir}/absolute_effect_diff_vs_maf_pb_chr{formatted_range}_py.png")
    plt.clf()

    # Plot a manhattan of absolute differences in betas
    merged_sorted = merged.sort_values(by=['CHR', 'POS'])
    merged_sorted['cumulative_position'] = merged_sorted.groupby('CHR').cumcount()
    # Create the Manhattan plot
    plt.figure(figsize=(15, 6))
    sns.scatterplot(x='cumulative_position', y='absolute_diff_betas', hue='CHR', data=merged_sorted, palette='viridis', alpha=0.7)
    # Add labels and title
    plt.xlabel('Cumulative Position')
    plt.ylabel('Absolute Difference Betas')
    plt.title('Manhattan Plot')
    plt.savefig(f"{outdir}/manhattan_abs_beta_exclusivity_chr{formatted_range}_py.png")
    print("Plotted Manhattan")

    # Looking at the top hit per gene (maybe we greater a greater hit rate per gene as opposed to variants?) [Make sure dfs are corrected same way]
    sc_res = sc_res.rename_axis(index={sc_res.index.names[0]: None})
    pb_res['qvalue_within_gene'] = pb_res.groupby('phenotype_id')['pval_nominal'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
    pb_res_min_qvalue_rows = pb_res.loc[pb_res.groupby('phenotype_id')['qvalue_within_gene'].idxmin()]
    pb_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(pb_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]
    sc_res['qvalue_within_gene'] = sc_res.groupby('Gene')['p.value'].transform(lambda x: smt.multipletests(x, method='fdr_bh')[1])
    sc_res_min_qvalue_rows = sc_res.loc[sc_res.groupby('Gene')['qvalue_within_gene'].idxmin()]
    sc_res_min_qvalue_rows['qvalue_across_genes'] = smt.multipletests(sc_res_min_qvalue_rows['qvalue_within_gene'], method='fdr_bh')[1]

    # Plotting the overlap of genes detected by each method at different resolutions
    pb_res_min_qvalue_rows.rename(columns={"phenotype_id": "Gene"}, inplace=True)
    merged_qvalue = sc_res_min_qvalue_rows.merge(pb_res_min_qvalue_rows, on="Gene", suffixes=('_saige', '_tensor'))
    merged_qvalue.to_csv(f"{tabdir}/min_qvalue_per_common_gene_saige_tensor_{formatted_range}.txt", sep = "\t", index=False) 
    cutoffs = [1, 0.05]
    for c in cutoffs:
        sub = merged_qvalue[(merged_qvalue['qvalue_across_genes_saige'] < c) | (merged_qvalue['qvalue_across_genes_tensor'] < c)]
        all_test = sub.shape[0]
        nsaige_win = sub[sub['qvalue_across_genes_saige'] < sub['qvalue_across_genes_tensor']].shape[0]
        perc_saige_win = "{:.2f}".format(100*(nsaige_win/all_test))
        ntensor_win = sub[sub['qvalue_across_genes_saige'] > sub['qvalue_across_genes_tensor']].shape[0]
        perc_tensor_win = "{:.2f}".format(100*(ntensor_win/all_test))
        print(f"At a p-value threshold of {c} across genes: Out of a possible {all_test} tests, SAIGE wins {nsaige_win} times ({perc_saige_win}%), Tensor wins {ntensor_win} times ({perc_tensor_win}%)")

    # Plot the effect across genes on a venn
    sub = merged_qvalue[(merged_qvalue['qvalue_across_genes_saige'] < 0.05) | (merged_qvalue['qvalue_across_genes_tensor'] < 0.05)]
    saige_hits = set(sub[sub['qvalue_across_genes_saige'] < 0.05]['Gene'].tolist())
    tensor_hits = set(sub[sub['qvalue_across_genes_tensor'] < 0.05]['Gene'].tolist())
    venn_labels = {'100': 'SAIGE-QTL', '010': 'TensorQTL', '110': 'Intersection'}

    plt.figure(figsize=(8, 6))
    venn2(subsets=[saige_hits, tensor_hits], set_labels=('SAIGE-QTL', 'TensorQTL'), set_colors=('skyblue', 'lightgreen'), alpha=0.7)
    plt.savefig(f"{outdir}/venn_saige_tensor_fdr_0.0_across_genes_chr{formatted_range}_py.png")
    plt.clf()

    # Also plot the correlation of betas/p-values for these significant effects - if only significant in one condition, plot one colour and another if the opposite
    sub = merged_qvalue[(merged_qvalue['qvalue_across_genes_saige'] < 0.05) | (merged_qvalue['qvalue_across_genes_tensor'] < 0.05)]
    max_effect = max([np.abs(max(sub['slope'].astype('float'))), np.abs(max(sub['BETA'].astype('float')))])
    sub['exclusivity'] = 'both'  # Default value for all rows
    sub.loc[(sub['qvalue_across_genes_saige'] < 0.05) & (sub['qvalue_across_genes_tensor'] > 0.05), 'exclusivity'] = 'SAIGE_only'
    sub.loc[(sub['qvalue_across_genes_saige'] > 0.05) & (sub['qvalue_across_genes_tensor'] < 0.05), 'exclusivity'] = 'Tensor_only'
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
    plt.text(max_effect, 0, f'Rho: {rho:.2f}\nP-value: {p_value:.4f}', color='blue', ha='right')
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

    # Plot the relationship between significance in one study and MAF
    sub['AF_Allele2'] = sub['AF_Allele2'].astype('float')
    sub['absolute_diff_betas'] = np.abs(sub['slope'] - sub['BETA'])
    sub['MAF'] = np.where(sub['AF_Allele2'] > 0.5, 1 - sub['AF_Allele2'], sub['AF_Allele2'])
    # Plot, colouring points by where they are significant
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='MAF', y='absolute_diff_betas', hue='exclusivity', data=sub)
    plt.xlabel('MAF')
    plt.ylabel('Absolute difference in effect size (SAIGE-QTL vs TensorQTL)')
    plt.xlim(0, 0.5)
    plt.title(f"FDR<0.05, Chromosome {formatted_range} only")
    plt.savefig(f"{outdir}/absolute_effect_diff_vs_maf_pb_fdr_0.05_chr{formatted_range}_py.png")
    plt.clf()

    # Plot the distribution of effect differences, within exclusivity
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    groups = np.unique(sub.exclusivity)
    for g in groups:
        sns.distplot(sub[sub['exclusivity'] == g].absolute_diff_betas, hist=False, rug=True, label=g)

    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.title('Distribution absolute differences in effect sizes / exclusivity')
    plt.xlabel('Difference in effect sizes')
    #plt.axvline(x = 0.7, color = 'red', linestyle = '--', alpha = 0.5)
    plt.savefig(f"{outdir}/dist_absolute_effect_diff_pb_fdr_0.05_chr{formatted_range}_py.png", bbox_inches='tight')
    plt.clf()

    # Plot the distribution of effect differences, within exclusivity
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    groups = np.unique(sub.exclusivity)
    for g in groups:
        sns.distplot(sub[sub['exclusivity'] == g].MAF, hist=False, rug=True, label=g)

    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.title('Distribution MAF / exclusivity')
    plt.xlabel('MAF')
    plt.xlim(0,0.5)
    plt.savefig(f"{outdir}/dist_maf_pb_fdr_0.05_chr{formatted_range}_py.png", bbox_inches='tight')
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

    # Plot examples for when SAIGE detects and Tensor does not - and vv
    adata = ad.read_h5ad(phenotype_file)
    adata = adata[adata.obs[aggregate_on] == level]

    # Load genotypes and subset for those with significant (across genes) effects in one or both
    saige_hit_only = merged_qvalue[(merged_qvalue['qvalue_across_genes_saige'] < 0.05) & (merged_qvalue['qvalue_across_genes_tensor'] > 0.05)].sort_values(by='qvalue_across_genes_saige', ascending=True)
    saige_hit_genes = saige_hit_only['Gene'].values
    tensor_hit_only = merged_qvalue[(merged_qvalue['qvalue_across_genes_tensor'] < 0.05) & (merged_qvalue['qvalue_across_genes_saige'] > 0.05)].sort_values(by='qvalue_across_genes_tensor', ascending=True).sort_values(by='qvalue_across_genes_tensor', ascending=True)
    tensor_hit_genes = tensor_hit_only['Gene'].values
    want_genes = np.concatenate([tensor_hit_genes, saige_hit_genes])
    want_vars = np.concatenate([tensor_hit_only['variant_id'].values, saige_hit_only['MarkerID'].values])
    G = read_plink1_bin(f"{geno_dir}/plink_genotypes_ordered.bed")
    G_id = pd.read_csv(f"{geno_dir}/plink_genotypes_ordered.bim", sep = "\t", header=None)
    G['variant'] = G_id.iloc[:,1]
    G = G[G.sample.isin(adata.obs.Corrected_genotyping_ID)]
    G = G[:,G['variant'].isin(want_vars)]

    # Plot the expression of these genes
    genotype_df = G.to_dataframe(name='genotype').reset_index()
    genotype_df = genotype_df.pivot(index='sample', columns='variant', values='genotype')
    expression_df = pd.DataFrame.sparse.from_spmatrix(adata.X)
    expression_df.columns=adata.var.index.values
    expression_df.index = adata.obs.index.values
    expression_df = expression_df.iloc[:,expression_df.columns.isin(want_genes)]
    genotype_df = genotype_df.reset_index()
    genotype_df = genotype_df.rename(columns={'sample': 'Corrected_genotyping_ID'})
    temp = adata.obs
    temp = temp.reset_index()
    genotype_df = genotype_df.merge(temp[["index", "Corrected_genotyping_ID"]], on="Corrected_genotyping_ID")
    genotype_df = genotype_df.set_index("index")
    expression_df = expression_df.merge(genotype_df, left_index=True, right_index=True)
    expression_df['sample'] = expression_df.index.str.split('-').str[2]
    expression_df['sample'] = expression_df['sample'].astype(str)

    # Plot some (saige hit only first)
    for r, row in saige_hit_only[:5].iterrows():
        try:
            marker = row['MarkerID']
            gene = row['Gene']
            print(f"Plotting association between {marker} and {gene}")
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
            ## Add text (NOTE: The pvalue plotted for the opposite test is not neccessarily from the same variant, but is the strongest assoc)
            saige_beta=row['BETA']
            saige_fdr=row['qvalue_across_genes_saige']
            tensor_beta=row['slope']
            tensor_fdr = row['qvalue_across_genes_tensor']
            text = f"SAIGE:beta = {float(saige_beta):.2f}, FDR={saige_fdr:.2e}. TENSOR (best assoc):{float(tensor_beta):.2f}, FDR={tensor_fdr:.2e}"
            ax.text(1, 1.1 * ax.get_ylim()[1], text, ha='center', va='bottom')
            # Add options 
            ax.set_xlabel(marker)
            ax.set_ylabel(gene)
            # save
            plt.savefig(f"{outdir}/saige_only_cis_{marker}_{gene}_sns.png")
            plt.clf()
        except Exception as e:
            # Print or log the error message
            print(f"Error processing row {r}: {e}")


    # Plot some (tensor hit only first)
    for r, row in tensor_hit_only[:5].iterrows():
        try:
            marker = row['variant_id']
            gene = row['Gene']
            print(f"Plotting association between {marker} and {gene}")
            ref = row['variant_id'].split(":")[-2]
            nonref = row['variant_id'].split(":")[-1]
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
            tensor_beta=row['slope']
            tensor_fdr=row['qvalue_across_genes_tensor']
            saige_beta=row['BETA']
            saige_fdr=row['qvalue_across_genes_saige']
            text = f"TENSOR:beta = {float(tensor_beta):.2f}, FDR={tensor_fdr:.2e}. SAIGE (best assoc):beta={float(saige_beta):.2f}, FDR={saige_fdr:.2e}"
            ax.text(1, 1.1 * ax.get_ylim()[1], text, ha='center', va='bottom')
            # Add options 
            ax.set_xlabel(marker)
            ax.set_ylabel(gene)
            # save
            plt.savefig(f"{outdir}/tensor_only_cis_{marker}_{gene}_sns.png")
            plt.clf()
        except Exception as e:
            # Print or log the error message
            print(f"Error processing row {r}: {e}")

    # Trans-analysis
    print("~~~~~~~~~~~~~~~~~~ Summarising the trans-eQTL analysis ~~~~~~~~~~~~~~~~~~")
    print("Loading in the data - all chromosomes")
    res = pd.read_csv(f"{catdir}/trans/all_nPC_{str(n_expr_pcs)}_trans_by_cis_no_cis_bonf.txt", sep = "\t")

    # Subset for the trans hits on genes tested by the pseudo-bulk analyses
    res = res[res['Gene'].isin(merged_qvalue['Gene'])]

    # Grab the top 5 and plot their trans-effects
    res_sig = res[res['bonf_qvalue_across_genes'] < 0.01]
    res_sig = res_sig.sort_values(by='p.value.bonf')
    res_plot = res_sig.head(50)

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
    circos.text(f"{level}\ntrans-eQTLs", deg=315, r=150, size=12)
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
        
    # Load the genotypes of these results
    G = read_plink1_bin(f"{geno_dir}/plink_genotypes_ordered.bed")
    G_id = pd.read_csv(f"{geno_dir}/plink_genotypes_ordered.bim", sep = "\t", header=None)
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

    for r, row in res_plot.iterrows():
        try:
            marker = row['MarkerID']
            gene = row['Gene']
            print(f"Plotting association between {marker} and {gene}")
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
            if len(plot_levels_ordered) != len(order):
                sns.violinplot(x=violin_df[marker], y=violin_df[gene], order=[2.0, 1.0], ax=ax)
            else:
                sns.violinplot(x=violin_df[marker], y=violin_df[gene], order=order, ax=ax)
            #
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
            plt.savefig(f"{outdir}/trans_{marker}_{gene}_trans_sns.png")
            plt.clf()
        #
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
    res['fdr_less_0.05'] = np.where(res['bonf_qvalue_across_genes'] < 0.05, 'Yes', 'No')
    res['MAF'] = np.where(res['AF_Allele2'].astype('float') > 0.5, 1 - res['AF_Allele2'].astype('float'), res['AF_Allele2'].astype('float'))
    res['logMean'] = np.log10(res['Mean'])

    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='MAF', y='logMean', data=res, hue="fdr_less_0.05")
    plt.xlabel('MAF')
    plt.ylabel('log10(Mean expression of gene)')
    plt.xlim(0, 0.5)
    plt.title(f"Detection of trans-eQTLs (MAF vs gene expression): FDR<0.05")
    plt.savefig(f"{outdir}/trans_detection_MAF_expr_py.png")
    plt.clf()

    # Plot a distribution of MAF across significance thresholds
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    groups = np.unique(res['fdr_less_0.05'])
    for g in groups:
        sns.distplot(res[res['fdr_less_0.05'] == g].MAF, hist=False, rug=True, label=g, kde_kws={'linewidth': 1.5})

    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.title('Distribution of trans-eQTL MAF across significance (FDR<0.01)')
    plt.xlabel('MAF')
    plt.xlim(0,0.5)
    #plt.axvline(x = 0.7, color = 'red', linestyle = '--', alpha = 0.5)
    plt.savefig(f"{outdir}/dist_trans_MAF_significance_py.png", bbox_inches='tight')
    plt.clf()

    # How does beta value vary with the MAF and gene expression? 
    res['absBETA'] = np.abs(res['BETA'].astype('float'))
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='MAF', y='absBETA', data=res, palette="viridis", alpha=0.8)
    plt.xlabel('MAF')
    plt.ylabel('absolute effect size (beta)')
    plt.xlim(0, 0.5)
    plt.title(f"Detection of trans-eQTLs (MAF vs absolute beta), colored by Mean expression")
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

    # Plot expression sparsity vs best beta/p-value
    column_names = [col for col in expression_df.columns if col.startswith('ENSG')]
    ensg_columns = expression_df[column_names]
    ensg_columns = ensg_columns.iloc[:,ensg_columns.columns.isin(res['Gene'])]
    is_zero = ensg_columns.eq(0)
    dense_zero_counts = is_zero.sparse.to_dense()
    gene_sparsity = pd.DataFrame(1-dense_zero_counts.sum()/dense_zero_counts.shape[0])
    gene_sparsity.reset_index(inplace=True)
    gene_sparsity.columns = ["Gene", "sparsity"]
    res = res.merge(gene_sparsity, how="left", on="Gene")

    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='sparsity', y='log10_p.value', hue="fdr_less_0.05", data=res)
    plt.xlabel('Gene sparsity (proportion non-zero counts)')
    plt.ylabel('log10(p-value)')
    plt.xlim(0, 1)
    plt.savefig(f"{outdir}/trans_detection_significance_sparsity_py.png")
    plt.clf()

    # Beta
    plt.figure(figsize=(8, 6))
    fig,ax = plt.subplots(figsize=(8,6))
    sns.scatterplot(x='sparsity', y='absBETA', hue="fdr_less_0.05", data=res)
    plt.xlabel('Gene sparsity (proportion non-zero counts)')
    plt.ylabel('absolute beta')
    plt.xlim(0, 1)
    plt.savefig(f"{outdir}/trans_detection_beta_sparsity_py.png")
    plt.clf()

if __name__ == '__main__':
    main()

    
    

