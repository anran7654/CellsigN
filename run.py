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
    """
    读取表达矩阵并做基础预处理
    """
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
            if cell_path is None or cell_path == "":
                cell_type = [s.rsplit("_", 1)[0] if "_" in s else s for s in cell_names]
                adata = sc.AnnData(
                    X=csr_matrix(data),
                    obs=pd.DataFrame({"celltype": cell_type}, index=cell_names),
                    var=pd.DataFrame(index=gene_names),
                )
            else:
                cell_dict = read_cell_type(cell_path)
                cell_type = []
                for cell in cell_names:
                    cell_type.append(cell_dict.get(cell))
                adata = sc.AnnData(
                    X=csr_matrix(data),
                    obs=pd.DataFrame({"celltype": cell_type}, index=cell_names),
                    var=pd.DataFrame(index=gene_names),
                )

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
    """
    备用的细胞聚类函数（当前脚本里不一定用得到）
    """
    sc.tl.pca(adata)
    data = adata.obsm["X_pca"]

    if num_cluster is None:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.louvain(adata, resolution=1.0)
        k = adata.obs["louvain"].nunique()
        print("Louvain cluster: " + str(k))

        if adata.n_obs < 2000:
            resolution = 0.8
        else:
            resolution = 0.5
        num_cluster = int(k * resolution) if int(k * resolution) >= 3 else 2

    agg_clustering = AgglomerativeClustering(n_clusters=num_cluster)
    adata.obs["celltype"] = agg_clustering.fit_predict(data)


def read_cell_type(filename):
    my_dict = {}
    with open(filename, "r") as file:
        for line in file:
            key, value = re.split("[\t,]", line.strip())
            my_dict[key] = value
    return my_dict


def get_top_genes(adata, cell_type, num_genes):
    gene_names = adata.uns["rank_genes_groups"]["names"][cell_type]
    log_fold_changes = adata.uns["rank_genes_groups"]["logfoldchanges"][cell_type]
    over_expressed_genes = [
        gene for gene, lfc in zip(gene_names, log_fold_changes) if lfc > 0
    ]
    top_gene = over_expressed_genes[:num_genes]
    return top_gene


def parse_args():
    parser = argparse.ArgumentParser(description="Main entrance of CellSigN")

    # 必需路径
    parser.add_argument("--file_path", type=str, required=True,
                        help="Path of gene expression data")
    parser.add_argument("--cell_path", type=str, required=True,
                        help="Path of cell labels")
    parser.add_argument("--pathway_file", type=str, required=True,
                        help="Path of pathway interaction file, e.g. Pathway_Commons.txt")
    parser.add_argument("--ligand_file", type=str, required=True,
                        help="Path of ligand-receptor interaction file")

    # 预处理与参数
    parser.add_argument("--min_cell", type=float, default=0.01,
                        help="Fraction of cells a gene must be expressed in (default: 0.01)")
    parser.add_argument("--min_gene", type=float, default=0.01,
                        help="Fraction of genes a cell must express (default: 0.01)")

    parser.add_argument("--no_normalize", action="store_true",
                        help="Disable library-size normalization")
    parser.add_argument("--no_log_trans", action="store_true",
                        help="Disable log1p transform")

    parser.add_argument("--cell_top_gene", type=int, default=500,
                        help="Top X genes for each celltype (default: 500)")
    parser.add_argument("--alpha", type=float, default=0.01,
                        help="Significant level, larger alpha leads to more edges (0-1, default: 0.01)")
    parser.add_argument("--max_length", type=int, default=10,
                        help="Maximum length of the receptor-TF pathway (default: 10)")

    parser.add_argument("--output_dir", type=str, default="Results",
                        help="Output directory for result files (default: Results)")
    parser.add_argument("--malignant_label", type=str, default="Malignant",
                        help="Celltype name for malignant cells in obs['celltype'] (default: Malignant)")

    return parser.parse_args()


def main():
    args = parse_args()

    print("Read single-cell data...")
    adata = read_file(
        file_path=args.file_path,
        cell_path=args.cell_path,
        min_cell=args.min_cell,
        min_gene=args.min_gene,
        normalize=not args.no_normalize,
        log_trans=not args.no_log_trans,
    )
    if adata is None:
        sys.exit(1)

    # 去掉没有 celltype 的细胞
    filtered_adata = adata[adata.obs["celltype"].notna() & (adata.obs["celltype"] != '""')].copy()
    adata = filtered_adata
    print("adata shape:", adata.shape)

    total_gene_dict = {var: idx for idx, var in enumerate(adata.var_names)}
    total_gene_list = adata.var_names
    data = adata.X.T.toarray()
    mean_all = np.mean(data, axis=1)
    variance_all = np.var(data, axis=1)

    print("Read ligand-receptor file...")
    lgrp_matrix = ligand_receptor(args.ligand_file, total_gene_dict, 0, 1, 2, 2)

    cell_types = adata.obs["celltype"].unique()
    print("Cell types:", cell_types)

    sc.tl.rank_genes_groups(
        adata,
        groupby="celltype",
        method="wilcoxon",
        n_genes=args.cell_top_gene,
    )

    dict_gene = {}
    tftg_final = {}
    pathway_final = {}
    interaction_type = defaultdict(set)

    # 输出预处理后的基因列表到 output_dir
    os.makedirs(args.output_dir, exist_ok=True)
    gene_list = list(adata.var_names)
    gene_list_file = os.path.join(args.output_dir, "processed_gene_list.txt")
    with open(gene_list_file, "w") as f:
        for gene in gene_list:
            f.write(gene + "\n")
    print(f"基因列表已保存到: {gene_list_file}")

    # 每个细胞类型跑一次 CSN + Pathway
    for cell_type in cell_types:
        cells_in_type = adata.obs["celltype"] == cell_type
        top_genes = get_top_genes(adata, cell_type, args.cell_top_gene)

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

        if sub_adata.n_obs > 2:
            gene_dict = {var: idx for idx, var in enumerate(sub_adata.var_names)}
            rp_tf, tf_tg, sub_type = pathway(
                args.pathway_file,
                gene_dict,
                col1_index=0,
                col2_index=1,
                col3_index=2,
            )
            for key, value in sub_type.items():
                interaction_type[key].update(value)

            rp_tf_adj, tf_tg_adj = csnet(
                sub_adata,
                mean_all2,
                variance_all2,
                adata.n_obs,
                cell_type,
                total_gene_dict,
                rt_matrix=rp_tf,
                tf_matrix=tf_tg,
                alpha=args.alpha,
            )

            pathway_final[cell_type] = rp_tf_adj
            tftg_final[cell_type] = tf_tg_adj

    malignant = args.malignant_label
    if malignant not in tftg_final or malignant not in pathway_final or malignant not in dict_gene:
        raise ValueError(
            f"Malignant label '{malignant}' not found in cell types: {list(cell_types)}"
        )

    tftg_malignant = tftg_final[malignant]
    pathway_malignant = pathway_final[malignant]
    malignant_gene = dict_gene[malignant]

    # 其它细胞 → Malignant
    for cell, gene_set in dict_gene.items():
        if cell == malignant:
            continue

        shortest_other_to_malignant = infer_pathway(
            gene_set, lgrp_matrix, pathway_malignant, tftg_malignant, args.max_length
        )
        if len(shortest_other_to_malignant) > 0:
            print(f"Found {len(shortest_other_to_malignant)} pathways from {cell} to {malignant}")
            output_file = os.path.join(args.output_dir, f"{cell}_to_{malignant}_pathway.txt")
            with open(output_file, "w") as f:
                f.write("Ligand\tReceptor\tMediator\tTF\tTarget\n")
                for sublist in shortest_other_to_malignant:
                    path_genes = [total_gene_list[i] for i in sublist]
                    col1 = path_genes[0]
                    col2 = path_genes[1]
                    col3 = ",".join(path_genes[2:-1])
                    col4 = path_genes[-1]
                    tgs = tftg_malignant.col[tftg_malignant.row == sublist[-1]].tolist()
                    col5 = ",".join([total_gene_list[i] for i in tgs])
                    f.write(f"{col1}\t{col2}\t{col3}\t{col4}\t{col5}\n")

        # Malignant → 其它细胞
        if cell in pathway_final and cell in tftg_final:
            adj_pathway = pathway_final[cell]
            adj_tftg = tftg_final[cell]
            shortest_malignant_to_other = infer_pathway(
                malignant_gene, lgrp_matrix, adj_pathway, adj_tftg, args.max_length
            )
            if len(shortest_malignant_to_other) > 0:
                print(f"Found {len(shortest_malignant_to_other)} pathways from {malignant} to {cell}")
                output_file = os.path.join(args.output_dir, f"{malignant}_to_{cell}_pathway.txt")
                with open(output_file, "w") as f:
                    f.write("Ligand\tReceptor\tMediator\tTF\tTarget\n")
                    for sublist in shortest_malignant_to_other:
                        path_genes = [total_gene_list[i] for i in sublist]
                        col1 = path_genes[0]
                        col2 = path_genes[1]
                        col3 = ",".join(path_genes[2:-1])
                        col4 = path_genes[-1]
                        tgs = adj_tftg.col[adj_tftg.row == sublist[-1]].tolist()
                        col5 = ",".join([total_gene_list[i] for i in tgs])
                        f.write(f"{col1}\t{col2}\t{col3}\t{col4}\t{col5}\n")

    print("CellSigN signaling inference finished.")


if __name__ == "__main__":
    main()
