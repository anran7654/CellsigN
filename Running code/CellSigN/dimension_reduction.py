import os
import sys
import argparse
import csv

import numpy as np
import scanpy as sc
import pandas as pd
import networkx as nx

from csn.cell_net import csnet
from data_process.read_expression import read_gene_expression
from data_process.read_interaction import pathway, ligand_receptor
from scipy.sparse import csr_matrix
from sklearn.cluster import AgglomerativeClustering
from csn.cell_communication import find_path_with_networkx, generate_concatenated_lists
from collections import defaultdict


def read_file(file_path,
              cell_path=None):
    if file_path.lower().endswith(".h5ad"):
        adata = sc.read(file_path)
    else:
        cell_names, gene_names, data = read_gene_expression(file_path)
        if cell_names is None:
            return None
        else:
            if len(cell_names) == data.shape[0] + 1:
                cell_names = cell_names[1:]
            elif len(cell_names) != data.shape[0]:
                print("Warning: cell count does not match expression matrix!")
                return None

            adata = sc.AnnData(X=csr_matrix(data),
                               obs=pd.DataFrame(index=cell_names),
                               var=pd.DataFrame(index=gene_names))

            cell_dict = read_cell_type(cell_path)

            adata_filtered = adata[adata.obs_names.isin(cell_dict.keys()), :]
            cell_names = adata_filtered.obs_names

            cell_type = []
            for cell in cell_names:
                cell_type.append(cell_dict.get(cell))

            adata_filtered.obs['celltype'] = cell_type

    if not adata_filtered.var_names.is_unique:
        adata_filtered.var_names_make_unique()

    return adata_filtered


def read_cell_type(filename):
    my_dict = {}
    with open(filename, 'r') as file:
        # 自动检测分隔符
        dialect = csv.Sniffer().sniff(file.read(1024))
        file.seek(0)  # 重置文件指针
        reader = csv.reader(file, delimiter=dialect.delimiter)
        for row in reader:
            if len(row) >= 2:
                key, value = row[0], row[1]
                my_dict[key] = value
    return my_dict


parser = argparse.ArgumentParser(description='Main entrance of DeepCSN')
parser.add_argument('--file_path', type=str, default='',
                    help='the path of gene expression data')
parser.add_argument('--cell_path', type=str, default='',
                    help='the path of cell labels')
args = parser.parse_args()

args.file_path = r'D:\Personal issues\CSNet\Expression data\CCCA\Data_Zhang2019_Ovarian/Data_Zhang2019_Ovarian.txt'
args.cell_path = r'D:\Personal issues\CSNet\Expression data\CCCA\Data_Zhang2019_Ovarian/Cells.csv'
n_clusters = 50

print('Read single-cell data...')
adata = read_file(args.file_path, args.cell_path)
if adata is None:
    sys.exit()

filtered_adata = adata[adata.obs['celltype'].notna() & (adata.obs['celltype'] != '""')].copy()
adata = filtered_adata
print(adata.shape)


cell_types = adata.obs['celltype'].unique()
print(cell_types)

list_data = []
for cell_type in cell_types:
    cell_indices = adata.obs[adata.obs['celltype'] == cell_type].index

    # If the number of cells for the current type is greater than 50
    if len(cell_indices) > 50:
        # Extract the gene expression matrix for the current cell type
        E = adata[cell_indices, :].X.toarray()

        # Apply AgglomerativeClustering to reduce the number of cells to 50
        clustering = AgglomerativeClustering(n_clusters=n_clusters)
        clustering_labels = clustering.fit_predict(E)  # Convert sparse matrix to dense

        # Compute the average expression for each cluster
        E_reduced = np.zeros((E.shape[1], n_clusters))
        for cluster_id in range(n_clusters):
            cluster_cells = np.where(clustering_labels == cluster_id)[0]
            avg_expr = E[cluster_cells].mean(axis=0)
            E_reduced[:, cluster_id] = avg_expr

        # Convert the reduced matrix into a pandas DataFrame
        reduced_df = pd.DataFrame(
            E_reduced,
            index=adata.var_names,
            columns=[f"{cell_type}_{i + 1}" for i in range(n_clusters)]
        )
        list_data.append(reduced_df)
        # Save the reduced matrix as a txt file
        # reduced_df.to_csv(f'D:/{cell_type}_expression.txt', sep='\t', header=True, index=True)
    else:
        E = adata[cell_indices, :].X.toarray().T
        reduced_df = pd.DataFrame(
            E,
            index=adata.var_names,
            columns=[f"{cell_type}_{i + 1}" for i in range(len(cell_indices))]
        )
        list_data.append(reduced_df)
        # reduced_df.to_csv(f'D:/{cell_type}_expression.txt', sep='\t', header=True, index=True)

merged_df = pd.concat(list_data, axis=1)
merged_df.to_csv(f'D:/Data_Zhang2019_Ovarian.txt', sep='\t', header=True, index=True)
