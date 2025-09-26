import os
import re
import sys
import argparse

import numpy as np
import scanpy as sc
import pandas as pd

from csn.cell_net import csnet
from data_process.read_expression import read_gene_expression
from data_process.read_interaction import pathway, ligand_receptor
from scipy.sparse import csr_matrix
from sklearn.cluster import AgglomerativeClustering
from csn.Infer_pathway import infer_pathway
from collections import defaultdict


def read_file(file_path,
              cell_path=None,
              min_cell=0.01,
              min_gene=0.01,
              normalize: bool = True,
              log_trans: bool = True):
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
            if cell_path is None:
                cell_type = [s.rsplit("_", 1)[0] if "_" in s else s for s in cell_names]
                adata = sc.AnnData(X=csr_matrix(data),
                                   obs=pd.DataFrame({'celltype': cell_type}, index=cell_names),
                                   var=pd.DataFrame(index=gene_names))
            else:
                cell_dict = read_cell_type(cell_path)
                cell_type = []
                for cell in cell_names:
                    cell_type.append(cell_dict.get(cell))
                adata = sc.AnnData(X=csr_matrix(data),
                                   obs=pd.DataFrame({'celltype': cell_type}, index=cell_names),
                                   var=pd.DataFrame(index=gene_names))
    if not adata.var_names.is_unique:
        adata.var_names_make_unique()

    # Filter cells
    if min_gene is not None:
        min_genes = max(1, int(min_gene * adata.n_vars))
        sc.pp.filter_cells(adata, min_genes=min_genes)

    # Filter genes
    if min_cell is not None:
        min_cells = max(1, int(min_cell * adata.n_obs))
        sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize and Log Transformation
    if normalize:
        sc.pp.normalize_total(adata)
    if log_trans:
        sc.pp.log1p(adata)

    return adata


def clustering_adata(adata, num_cluster: int = None):
    sc.tl.pca(adata)
    data = adata.obsm['X_pca']

    if num_cluster is None:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # use pcs to compute the neighborhood graph
        sc.tl.louvain(adata, resolution=1.0)
        k = adata.obs['louvain'].nunique()
        print('Louvain cluster: ' + str(k))

        if adata.n_obs < 2000:
            resolution = 0.8
        else:
            resolution = 0.5
        num_cluster = int(k * resolution) if int(k * resolution) >= 3 else 2

    agg_clustering = AgglomerativeClustering(n_cluster=num_cluster)
    adata.obs['celltype'] = agg_clustering.fit_predict(data)


def read_cell_type(filename):
    my_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            key, value = re.split('[\t,]', line.strip())
            my_dict[key] = value
    return my_dict


def get_top_genes(adata, cell_type, num_genes):
    gene_names = adata.uns['rank_genes_groups']['names'][cell_type]
    log_fold_changes = adata.uns['rank_genes_groups']['logfoldchanges'][cell_type]
    over_expressed_genes = [
        gene for gene, lfc in zip(gene_names, log_fold_changes) if lfc > 0
    ]
    top_gene = over_expressed_genes[:num_genes]
    return top_gene


parser = argparse.ArgumentParser(description='Main entrance of CellSign')
parser.add_argument('--file_path', type=str, default='',
                    help='the path of gene expression data')
parser.add_argument('--cell_path', type=str, default='',
                    help='the path of cell labels')
parser.add_argument('--min_cell', type=float, default=0.01,
                    help='parameter for gene filtering')
parser.add_argument('--min_gene', type=float, default=0.01,
                    help='parameter for cell filtering')
parser.add_argument('--normalize', type=bool, default=True,
                    help='normalize cells')
parser.add_argument('--log_trans', type=bool, default=True,
                    help='logarithm expression')
parser.add_argument('--cell_top_gene', type=int, default=500,
                    help='take the top X genes expressed in each cell (default: 1000)')
parser.add_argument('--alpha', type=float, default=0.01,
                    help='Significant level, larger alpha leads to more edges (default: 0.01), 0-1')
parser.add_argument('--max_length', type=int, default=10,
                    help='Maximum length of the receptor-TF pathway')
args = parser.parse_args()

args.file_path = r'D:\重生之我是研究生\代码\Evaluation_data\Data_Wang2019_Brain_SF4297\Data_Wang2019_Brain_SF4297.txt'
args.cell_path = r'D:\重生之我是研究生\代码\Evaluation_data\Data_Wang2019_Brain_SF4297/cell_type.csv'
pathway_file = r'D:\重生之我是研究生\代码\interaction data\human/Pathway_Commons.txt'
ligand_file = r'D:\重生之我是研究生\代码\interaction data\human/ligand-receptor.txt'

print('Read single-cell data...')
adata = read_file(args.file_path, args.cell_path, min_cell=None, min_gene=None,
                  normalize=args.normalize, log_trans=args.log_trans)
if adata is None:
    sys.exit()

filtered_adata = adata[adata.obs['celltype'].notna() & (adata.obs['celltype'] != '""')].copy()
adata = filtered_adata
print(adata.shape)


# sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg_top_genes)
# hvg_genes = adata.var[adata.var['highly_variable']].index
# adata = adata[:, list(hvg_genes)]

total_gene_dict = {var: idx for idx, var in enumerate(adata.var_names)}
total_gene_list = adata.var_names
data = adata.X.T.toarray()
mean_all = np.mean(data, axis=1)  # Calculate the mean of each row
variance_all = np.var(data, axis=1)  # Calculate the variance of each row

print('Read ligand-receptor file...')
lgrp_matrix = ligand_receptor(ligand_file, total_gene_dict, 0, 1, 2, 2)

cell_types = adata.obs['celltype'].unique()
print(cell_types)
sc.tl.rank_genes_groups(adata, groupby='celltype', method='wilcoxon', n_genes=args.cell_top_gene)

dict_gene = {}
tftg_final = {}
pathway_final = {}
interaction_type = defaultdict(set)  # Use defaultdict to simplify merging

# 将预处理后的基因列表保存到指定文件
gene_list = list(adata.var_names)  # 提取基因列表
output_file = r"D:\重生之我是研究生\代码\Evaluation_data\Data_Wang2019_Brain_SF4297/processed_gene_list.txt"  # 输出文件路径

with open(output_file, 'w') as f:
    for gene in gene_list:
        f.write(gene + '\n')  # 写入文件
print(f"基因列表已保存到: {output_file}")


for cell_type in cell_types:
    cells_in_type = adata.obs['celltype'] == cell_type
    top_genes = get_top_genes(adata, cell_type, args.cell_top_gene)  # str格式的list，存储基因名
    mean_all2 = {}
    variance_all2 = {}
    gene_set = set()
    for key in top_genes:
        gene_index = total_gene_dict[key]
        gene_set.add(gene_index)
        mean_all2[key] = mean_all[gene_index]
        variance_all2[key] = variance_all[gene_index]
    dict_gene[cell_type] = gene_set

    sub_adata = adata[cells_in_type, adata.var_names.isin(top_genes)].copy()
    # print(sub_adata)

    if sub_adata.n_obs > 2:
        gene_dict = {var: idx for idx, var in enumerate(sub_adata.var_names)}
        rp_tf, tf_tg, sub_type = pathway(pathway_file, gene_dict, col1_index=0, col2_index=1, col3_index=2)
        for key, value in sub_type.items():
            interaction_type[key].update(value)  # Merge sets with update()

        rp_tf_adj, tf_tg_adj = csnet(sub_adata, mean_all2, variance_all2, adata.n_obs, cell_type,
                                     total_gene_dict, rt_matrix=rp_tf, tf_matrix=tf_tg, alpha=args.alpha)
        pathway_final[cell_type] = rp_tf_adj
        tftg_final[cell_type] = tf_tg_adj

os.makedirs('Results', exist_ok=True)
file_path = f"Results/"
tftg_malignant = tftg_final['Malignant']  # 跟scATAC-seq数据的邻接矩阵取交集
pathway_malignant = pathway_final['Malignant']
malignant_gene = dict_gene['Malignant']
for cell, gene_set in dict_gene.items():
    if cell != 'Malignant':
        shortest_other_to_malignant = infer_pathway(gene_set, lgrp_matrix, pathway_malignant, tftg_malignant, args.max_length)
        if len(shortest_other_to_malignant) > 0:
            print(f"Found {len(shortest_other_to_malignant)} pathways from {cell} to Malignant")
            output_file = f'Results/{cell}_to_Malignant_pathway.txt'
            with open(output_file, 'w') as f:
                f.write(f"Ligand\tReceptor\tMediator\tTF\tTarget\n")
                for sublist in shortest_other_to_malignant:
                    pathway = [total_gene_list[i] for i in sublist]
                    col1 = pathway[0]  # First element
                    col2 = pathway[1]  # Second element
                    col3 = ",".join(pathway[2:-1])  # Elements between second and last
                    col4 = pathway[-1]  # Last element
                    tgs = tftg_malignant.col[tftg_malignant.row == sublist[-1]].tolist()
                    col5 = ",".join([total_gene_list[i] for i in tgs])
                    f.write(f"{col1}\t{col2}\t{col3}\t{col4}\t{col5}\n")

        if cell in pathway_final and cell in tftg_final:
            adj_pathway = pathway_final[cell]
            adj_tftg = tftg_final[cell]
            shortest_malignant_to_other = infer_pathway(malignant_gene, lgrp_matrix, adj_pathway, adj_tftg,
                                                        args.max_length)
            if len(shortest_malignant_to_other) > 0:
                print(f"Found {len(shortest_malignant_to_other)} pathways from Malignant to {cell}")
                output_file = f'Results/Malignant_to_{cell}_pathway.txt'
                with open(output_file, 'w') as f:
                    f.write(f"Ligand\tReceptor\tMediator\tTF\tTarget\n")
                    for sublist in shortest_malignant_to_other:
                        pathway = [total_gene_list[i] for i in sublist]
                        col1 = pathway[0]  # First element
                        col2 = pathway[1]  # Second element
                        col3 = ",".join(pathway[2:-1])  # Elements between second and last
                        col4 = pathway[-1]  # Last element
                        tgs = adj_tftg.col[adj_tftg.row == sublist[-1]].tolist()
                        col5 = ",".join([total_gene_list[i] for i in tgs])
                        f.write(f"{col1}\t{col2}\t{col3}\t{col4}\t{col5}\n")

'''
total_rptf = {}
final_paths = {}
for type1 in cell_types:
    for type2 in cell_types:
        if type1 != type2 and type1 in total_tftg and type2 in total_tftg:
            tf_tg1 = total_tftg[type1]
            tg_tf1 = total_tgtf[type1]
            tf_tg2 = total_tftg[type2]
            rp_tf_adj2 = total_rptf_adj[type2]

            rp_tf = {}
            runned = set()
            for tg1 in tg_tf1.keys():
                if tg1 in total_lgrp:
                    rp_set = total_lgrp[tg1]
                    for rp in rp_set:
                        if rp not in runned:
                            runned.add(rp)
                            if np.any(rp_tf_adj2[rp] != 0):
                                rptf_list = []
                                for tf in tf_tg2.keys():
                                    # path = bfs_shortest_path(rp_tf_adj2, rp, tf, max_length=args.max_length)
                                    path = find_path_with_networkx(rp_tf_adj2, rp, tf, max_length=args.max_length)
                                    if path is not None:
                                        rptf_list.append(path)
                                rp_tf[rp] = rptf_list

            if len(rp_tf) > 0:
                result = generate_concatenated_lists(tf_tg1, total_lgrp, rp_tf, tf_tg2)
                if len(result) > 0:
                    os.makedirs('Results', exist_ok=True)
                    file_path = f"Results/{type1} to {type2} pathways.txt"
                    with open(file_path, 'w') as f:
                        f.write(f"{type1}_TF\t{type1}_target\t{type2}_pathway\t{type2}_target\n")
                        for sublist in result:
                            # Process each element in the sublist
                            line_elements = []
                            for element in sublist:
                                if isinstance(element, int):
                                    line_elements.append(total_gene_list[element])
                                elif isinstance(element, list):
                                    line = []
                                    for name in element:
                                        line.append(total_gene_list[name])
                                    line_elements.append(','.join(map(str, line)))
                                elif isinstance(element, set):
                                    line = []
                                    for name in element:
                                        line.append(total_gene_list[name])
                                    line_elements.append(','.join(map(str, line)))
                                else:
                                    raise ValueError(f"Unsupported element type: {type(element)}")

                            # Join all elements of the line with tabs and write to file
                            f.write('\t'.join(map(str, line_elements)) + '\n')
                    print(f"Signal transduction from {type1} to {type2} complete")
'''