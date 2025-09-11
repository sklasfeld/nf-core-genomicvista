#!/usr/bin/env python3

import argparse
from tqdm import tqdm
import json
import pandas as pd
import numpy as np
import h5py
from pathlib import Path
from scipy.sparse import csr_matrix
import anndata as ad
import umap
import zarr

def count_lines(filepath):
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f)
    
def count_columns(filepath, sep='\t', skiprows=0):
    with open(filepath, 'r') as f:
        for _ in range(skiprows):
            next(f)
        first_line = next(f)
    return len(first_line.strip().split(sep))

def parse_gtf_attributes(attr_string):
    attributes = {}
    for item in attr_string.strip().split(';'):
        if item.strip():
            parts = item.strip().split(' ', 1) # Split only on the first space
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().strip('"') # Remove quotes
                attributes[key] = value
    return attributes

# function to generate subject ID from sampID
def getSubID(sampID):
    samp_arr = sampID.split('-')
    subID = '-'.join(samp_arr[0:2])
    return (subID)

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess expression data')
    parser.add_argument('--expression-matrix', required=True, help='Path to expression matrix')
    parser.add_argument('--tissue-attributes', required=True, help='Path to tissue attributes')
    parser.add_argument('--subject-attributes', required=True, help='Path to subject attributes')
    parser.add_argument('--annotation-gtf', required=True, help='Path to annotation GTF file')
    parser.add_argument('--output', required=True, help='Output processed data file')
    parser.add_argument('--stats-output', required=True, help='Output statistics file')
    parser.add_argument('--chunk_size', required=False, help='Chunk size for processing `expression-matrix`', default=1000, type=int)
    parser.add_argument('--skip_expression_rows', required=False, help='If the `expression-matrix` has metadata above the header, set this to the number of lines of meta data to skip at start of `expression-matrix`. Default is 2.', default=2, type=int)
    parser.add_argument('--expression_meta_cols', required=False, nargs='*', default=2, help='Number of columns in the `expression-matrix` file with feature information such as gene or transcript name. Default is 2 (Name & Description).', default=2, type=int)
    return parser.parse_args()

def load_data(expression_file, tissue_file, subject_file, annotation_gtf, chunk_size=1000, skip_expression_rows=2, expression_meta_cols=2):
    
    """Load all input data files"""

    print(f"Loading expression data from {expression_file}")
    # get the matrix dimensions
    nrows = count_lines(expression_file) - skip_expression_rows
    ncols = count_columns(expression_file, skiprows=skip_expression_rows)
    total_iterations = np.ceil(nrows / chunk_size)

    # first get the gene/transcript meta data from the expression data file
    if expression_meta_cols > 1:
        genes_meta_df = pd.concat(
            [chunk for chunk in tqdm(pd.read_csv(expression_file, sep="\t", skiprows=skip_expression_rows, usecols=list(range(0,expression_meta_cols)), compression='gzip', chunk_size=chunk_size), desc='Loading gene data', total=total_iterations)]
        )
        genes_meta_df.set_index(genes_meta_df.columns[0], inplace=True)
    else:
        genes_meta_df = pd.DataFrame()
    
    
    # Then get the TPMs from the expression data file
    expression_data = pd.concat(
        [chunk for chunk in tqdm(pd.read_csv(expression_file, sep="\t", skiprows=skip_expression_rows, usecols=range(expression_meta_cols,expression_meta_cols+ncols), compression='gzip', dtype=np.float64, chunk_size=chunk_size), desc='Loading expression data', total=total_iterations)]
    )

    print(f"Loading gene annotations from {annotation_gtf}")
    gtf_df = pd.read_csv(annotation_gtf, sep='\t', comment='#', header=None,dtype=str)
    gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gene_gtf_df = gtf_df[gtf_df["feature"] == "gene"].copy()
    add_columns = ['gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_type', 'transcript_name']
    gene_gtf_df.loc[:, add_columns] = None
    for index, row in gene_gtf_df.iterrows():
        attributes_dict = parse_gtf_attributes(row['attribute'])
    #gene_gtf_df.at[index,'attribute'] = attributes_dict
    for col in add_columns:
        if col in attributes_dict.keys():
            gene_gtf_df.loc[index,col] = attributes_dict[col]
   
    gene_gtf_df = gene_gtf_df.set_index('gene_id').copy()
    gene_gtf_df.index.name =None
    if len(genes_meta_df) > 0:
        genes_meta_df.index.name =None
        attribute_info = genes_meta_df.merge(gene_gtf_df, left_index=True, right_index=True, how='left')
    else:
        attribute_info = gene_gtf_df

    print(f"Loading tissue attributes from {tissue_file}")

    nrows = count_lines(tissue_file)
    total_iterations = np.ceil(nrows / chunk_size)
    tissue_attrs = pd.concat(
        [chunk for chunk in tqdm(pd.read_csv(tissue_file, sep="\t", chunk_size=chunk_size), desc='Loading subject data', total=total_iterations)]
    )
    # FOR NOW, WE WILL ONLY KEEP INFORMATION ABOUT THE TISSUES SAMPLED (in the future we can include more columns if needed)
    tissue_attrs = tissue_attrs.loc[tissue_attrs["SMAFRZE"]=="RNASEQ",["SAMPID", "SMTS", "SMTSD"]].copy()
    tissue_attrs = (
        tissue_attrs.rename(
            columns={
                "SMTS": "TISSUE",
                "SMTSD": "TISSUE_DETAILED"
            }
        )
    )

    # Add Subject IDs to tissue attributes
    tissue_attrs["SUBJID"] = tissue_attrs["SAMPID"].apply(getSubID)
    # Set SAMPID as index
    tissue_attrs.set_index("SAMPID", inplace=True)
    tissue_attrs = tissue_attrs.loc[expression_data.columns.tolist(),:]
    
    print(f"Loading subject attributes from {subject_file}")
    subject_attrs = pd.read_csv(subject_file, sep="\t", index_col=0)

    # Build AnnData object
    print("Building AnnData object")
    counts = csr_matrix(
        expression_data.transpose(), 
        dtype=np.float32)
    adata = ad.AnnData(counts)

    # Now, we provide the index to both the `var` (transcript/gene info) and `obs` (tissue/subject info) axes 
    adata.var_names = attribute_info.index.tolist()
    adata.var = attribute_info
    adata.obs_names = expression_data.columns.tolist()
    adata.obs["TISSUE"] = pd.Categorical(tissue_attrs['TISSUE'])  
    adata.obs["TISSUE_DETAILED"] = pd.Categorical(tissue_attrs['TISSUE_DETAILED']) 
    adata.obs["SUBJID"] = tissue_attrs['SUBJID']
    adata.obs = adata.obs.join(subject_attrs, on="SUBJID", how="left")
    
    return adata

def preprocess_data(adata):
    """Preprocessing steps"""
    
    

    # TODO: Add your preprocessing steps here
    # Examples:
    # - Quality control filtering
    # - Normalization
    # - Outlier detection
    # - Data transformation
    
    # log transformation of the expression data
    adata.layers["log_transformed"] = np.log1p(adata.X)

    # generate a UMAP embedding of the log1p transformed data
    reducer = umap.UMAP(random_state=42)
    adata.obsm["Xlog_umap"] = reducer.fit_transform(adata.layers["log_transformed"])

    return (adata)

def save_processed_data(adata, output_file):
    """Save processed data to HDF5 format"""
    adata.write_zarr(output_file)
    

def main():
    args = parse_args()
    
    # Load data into AnnData object
    adata = load_data(
        args.expression_matrix, 
        args.tissue_attributes, 
        args.subject_attributes,
        args.annotation_gtf,
        chunk_size=args.chunk_size,
        skip_expression_rows=args.skip_expression_rows,
        expression_meta_cols=args.expression_meta_cols
    )
    
    # Preprocess the data
    adata = preprocess_data(adata)
    
    # Save results
    save_processed_data(adata, args.output)
    
    # Save summary
    with open(args.summary_output, 'w') as f:
        f.write("Data Preprocessing Summary\n")
        f.write("=" * 30 + "\n")
        f.write("\nAnnData Object Structure\n")
        f.write("=" * 25 + "\n")
        f.write(f"Observations (samples): {adata.n_obs}\n")
        f.write(f"Variables (genes): {adata.n_vars}\n")
        f.write(f"Observation metadata columns: {list(adata.obs.columns)}\n")
        f.write(f"Variable metadata columns: {list(adata.var.columns)}\n")
        
        if adata.raw is not None:
            f.write(f"Raw data preserved: {adata.raw.n_vars} genes\n")
    
    print("Preprocessing completed successfully!")
    print(f"AnnData object: {adata.n_obs} samples × {adata.n_vars} genes")
    print(f"Summary saved to: {args.summary_output}")

if __name__ == '__main__':
    main()