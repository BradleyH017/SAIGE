#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-23'
__version__ = '0.0.1'

####### Python script to prepare input files for SAIGE

# Load libraries
##################################################
#### Bradley August 2023
# Atlassing the rectum scRNAseq
# Using all of the expression data, not just limited to that which is genotyped.
##################################################

# Load packages
import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mp
import kneed as kd
import math
import scipy.stats as st
from scipy import stats
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy.optimize import fsolve
from scipy.optimize import brentq
import argparse
import kneed as kd
from numpy import asarray
from numpy import savetxt
print("Loaded libraries")

# Define covariate process function
def preprocess_covariates(df, scale_covariates):
    processed_df = df.copy()
    for column in df.columns:
        if(scale_covariates == "true"):
            if df[column].dtype == 'float64':  # Check if the column is of type float
                # Scale and center numerical columns
                scaler = StandardScaler()
                processed_df[column] = scaler.fit_transform(df[[column]])
        elif pd.api.types.is_categorical_dtype(df[column]):  # Check if the column is categorical
            # Create dummy variables for categorical columns
            dummy_columns = pd.get_dummies(df[column], prefix=column, drop_first=True)
            # Remove the original categorical column
            processed_df.drop(column, axis=1, inplace=True)
            # Replace the dummy variable column name with the original column name
            dummy_columns.columns = [column]
            # Concatenate the dummy columns to the processed DataFrame
            processed_df = pd.concat([processed_df, dummy_columns], axis=1)
    return processed_df

# Inherit option
def parse_options():    
    # Inherit options
    parser = argparse.ArgumentParser(
            description="""
                Prepping files for the SAIGEQTL analysis
                """
        )

    parser.add_argument(
            '-ph', '--phenotype__file',
            action='store',
            dest='phenotype__file',
            required=True,
            help=''
        )

    parser.add_argument(
            '-a', '--aggregate_on',
            action='store',
            dest='aggregate_on',
            required=True,
            help=''
        )

    parser.add_argument(
            '-g', '--genotype_pc__file',
            action='store',
            dest='genotype_pc__file',
            required=True,
            help=''
        )

    parser.add_argument(
            '-id', '--genotype_id',
            action='store',
            dest='genotype_id',
            required=True,
            help=''
        )

    parser.add_argument(
            '-s', '--sample_id',
            action='store',
            dest='sample_id',
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
            '-p', '--nperc',
            action='store',
            dest='nperc',
            required=True,
            help=''
        )

    parser.add_argument(
            '-col', '--condition_col',
            action='store',
            dest='condition_col',
            required=False,
            default='',
            help=''
        )

    parser.add_argument(
            '-cond', '--condition',
            action='store',
            dest='condition',
            required=False,
            default='',
            help=''
        )

    parser.add_argument(
            '-covs', '--covariates',
            action='store',
            dest='covariates',
            required=True,
            help=''
        )

    parser.add_argument(
            '-xpca', '--expression_pca',
            action='store',
            dest='expression_pca',
            required=True,
            help=''
        )

    parser.add_argument(
            '-sc', '--scale_covariates',
            action='store',
            dest='scale_covariates',
            required=True,
            help=''
        )
    
    parser.add_argument(
            '-l', '--level',
            action='store',
            dest='level',
            required=True,
            help=''
        )
    
    return parser.parse_args()



# Define the main script
def main():
    inherited_options = parse_options()
    phenotype__file = inherited_options.phenotype__file
    aggregate_on = inherited_options.aggregate_on
    genotype_pc_file = inherited_options.genotype_pc__file
    genotype_id = inherited_options.genotype_id
    sample_id = inherited_options.sample_id
    general_file_dir = inherited_options.general_file_dir
    nperc = int(inherited_options.nperc)
    condition_col = inherited_options.condition_col
    condition = inherited_options.condition
    covariates = inherited_options.covariates
    expression_pca = inherited_options.expression_pca
    scale_covariates = inherited_options.scale_covariates
    level = inherited_options.level
    
    print(f"~~~~~~~~~~~~Working on: {level}~~~~~~~~~~~~~~~~")

    # Load in the adata object
    print("Loading object")
    adata = ad.read_h5ad(phenotype__file)

    # Load the genotype PCs
    print("Loading genotype PCs")
    geno_pcs = pd.read_csv(genotype_pc_file, sep = "\t")
    geno_pcs = geno_pcs.set_index("IID")
    geno_pcs = geno_pcs.iloc[:,1:]
    geno_pcs.reset_index(inplace=True)
    geno_pcs.rename(columns={"IID": genotype_id}, inplace=True)

    # Subset for the cells we want here (based on the condition column and value specified)
    if condition_col != "NULL":
        print("Subsetting for the condition")
        adata = adata[adata.obs[condition_col] == condition]

    # Replace the ages of those missing [SPECIFIC TO OUR DATA/SAMPLES]
    adata.obs.loc[adata.obs['sanger_sample_id'].isin(['OTARscRNA9294497', 'OTARscRNA9294498']), 'age_imputed'] = 56.5
    covs_use=covariates.split(',')

    # Define savedir
    if condition_col != "NULL":
        savedir=f"{general_file_dir}/{aggregate_on}/{level}/{condition_col}/{condition}"
    else:
        savedir=f"{general_file_dir}/{aggregate_on}/{level}"
    
    print(f"Files will be saved in {savedir}")
    
    
    if os.path.exists(savedir) == False:
        os.makedirs(savedir, exist_ok=True)
    print("Filtering anndata")
    temp = adata[adata.obs[aggregate_on] == level]
    
    
    # Filter for intersection with genotypes
    temp = temp[temp.obs[genotype_id].isin(geno_pcs[genotype_id])]
    
    
    # Identify genes expressed in > n% of cells
    print("Filtering lowly expressed genes")
    counts=temp.X
    count_per_column = np.array((counts > 1).sum(axis=0))[0]
    min_cells = math.ceil(temp.shape[0]*(int(nperc)/100))
    keep_genes = np.where(count_per_column > counts.shape[0] * nperc/100)[0]
    keep_genes_names = temp.var.index[keep_genes].values
    # Subset counts for these genes
    counts = counts[:,keep_genes]
    print(f"Final shape is:{counts.shape}")
    
    
    # Preprocess the covariates (scale continuous, dummy for categorical)
    print("Extracting and sorting covariates")
    to_add = temp.obs[covs_use]
    # Preprocess the covariates (scale continuous, dummy for categorical)
    to_add = preprocess_covariates(to_add, scale_covariates)
    # Bind this onto the counts
    counts = pd.DataFrame.sparse.from_spmatrix(counts)
    counts.index = temp.obs.index
    counts.columns = keep_genes_names
    counts = counts.merge(to_add, left_index=True, right_index=True)
    # Add the donor ID (genotyping ID so that we match the genotypes)
    counts[genotype_id] = temp.obs[genotype_id]
    # Also add the genotyping PCs
    index_use = counts.index
    counts = counts.merge(geno_pcs, on=genotype_id, how='left')
    counts.set_index(index_use, inplace=True)
    
    
    # If a covariate has ':' in it's name, this will throw errors in SAIGE, replace this with '_'
    counts.columns=counts.columns.str.replace(':', '_')
    
    
    # Compute and add the expression PCs
    if expression_pca == "true":
        print("Computing expression PCs")
        sc.pp.normalize_total(temp, target_sum=1e4)
        sc.pp.log1p(temp)
        sc.pp.highly_variable_genes(temp, flavor="seurat", n_top_genes=2000)
        sc.pp.scale(temp, max_value=10)
        sc.tl.pca(temp, svd_solver='arpack')
        pca_variance=pd.DataFrame({'x':list(range(1, 51, 1)), 'y':temp.uns['pca']['variance']})
        # Identify 'knee'
        knee=kd.KneeLocator(x=list(range(1, 51, 1)), y=temp.uns['pca']['variance'], curve="convex", direction = "decreasing")
        knee_point = knee.knee
        print('Knee: ', knee_point)
        # Save knee
        np.savetxt(f"{savedir}/knee.txt", [knee_point], delimiter=',', fmt='%s')
        # Append the PC matrix onto the count data (up to 20)
        loadings = pd.DataFrame(temp.obsm['X_pca'])
        loadings = loadings.iloc[:,0:20]
        loadings.index = temp.obs.index
        loadings.rename(columns=lambda x: f'xPC{x+1}', inplace=True)
        counts = counts.merge(loadings, left_index=True, right_index=True)
        
        
    # Save this as a dataframe in a directory specific to this resolution - make this if not already
    print("Saving")
    # Save list of gene names to test
    with open(f"{savedir}/test_genes.txt", 'w') as file:
        for item in keep_genes_names:
            file.write(f"{item}\n")
            
            
    # Save counts + covariates
    # NOTE: This currently isn't implemented to save sparse or compressed files (as SAIGE requires dense input) - 80k cells and 10k genes = ~2.5GB file
    counts.to_csv(f"{savedir}/saige_filt_expr_input.txt", sep = "\t", index=False)
    

if __name__ == '__main__':
    main()
    
