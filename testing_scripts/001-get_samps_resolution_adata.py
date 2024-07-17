#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-27'
__version__ = '0.0.1'

# Load packages
import anndata as ad
import scanpy as sc
import argparse
import numpy as np
import os

# Inherit options
def parse_options():  
    parser = argparse.ArgumentParser(
            description="""
                Script to perform scANVI integration of 
                """
        )

    parser.add_argument(
            '-p', '--phenotype__file',
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
            '-o', '--general_file_dir',
            action='store',
            dest='general_file_dir',
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
    
    return parser.parse_args()


# Define the main script
def main():
    inherited_options = parse_options()
    phenotype__file = inherited_options.phenotype__file
    aggregate_on = inherited_options.aggregate_on
    genotype_id = inherited_options.genotype_id
    general_file_dir = inherited_options.general_file_dir
    
    # Load the phenotype file
    adata=ad.read_h5ad(phenotype__file)

    # Extract the category levels
    cats = np.unique(adata.obs[aggregate_on])
    cats

    # Save this in the general file dir, making this directory uf not already present
    general_file_dir_cat=f"{general_file_dir}/{aggregate_on}"
    if os.path.exists(general_file_dir_cat) == False:
        os.makedirs(general_file_dir_cat, exist_ok=True)


    with open(f"{general_file_dir_cat}/levels.txt", 'w') as file:
            for item in cats:
                file.write(f"{item}\n")


    # Save a list of samples
    # Insert any code here if samples have missing covariates etc
    with open(f"{general_file_dir}/samples.txt", 'w') as file:
            for item in np.unique(adata.obs[genotype_id]):
                file.write(f"{item}\n")

if __name__ == '__main__':
    main()
# DONE!