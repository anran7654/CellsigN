import os
import numpy as np
from scipy.stats import norm, f, pearsonr, spearmanr, kendalltau
from scipy.sparse import csr_matrix, coo_matrix
from csn.correlation_test import pearson_correlation, statistic, var_net
from sklearn.decomposition import PCA
from collections import defaultdict


def pearson_net(data, p_value):
    n1 = data.shape[1]
    '''
    if n1 > 50:
        pca = PCA(n_components=50)
        data = pca.fit_transform(data)
        n1 = data.shape[1]  # n1是列数
    '''
    pearson_matrix = pearson_correlation(data)
    csn1 = statistic(pearson_matrix, n1, p_value)
    # spearman_matrix = spearman_correlation(data)
    # csn2 = statistic(spearman_matrix, n1, p_value)
    # union = (csn1 + csn2).astype(bool).astype(int)
    return csr_matrix(csn1)


def variance_direct(dif_array, gene_list, mean_all, variance_all, n1, mean_cell, variance_cell, n2, p_value):
    n3 = len(variance_all)
    sig_gene = var_net(list(variance_all.values()), n1 - 1, variance_cell, n2 - 1, p_value)
    csn = np.zeros((n3, n3)).astype(bool)

    rows, cols = dif_array.nonzero()
    for i, j in zip(rows, cols):
        gene1 = gene_list[i]
        gene2 = gene_list[j]
        if sig_gene[i] > 0 and sig_gene[j] > 0:
            combined_variance1 = (((n1 - 1) * variance_all[gene1] + (n1 - 1) * variance_all[gene2] + (n1 * n1) / (n1 + n1)
                                   * (mean_all[gene1] - mean_all[gene2]) ** 2) / (n1 + n1 - 1))
            combined_variance2 = (((n2 - 1) * variance_cell[i] + (n2 - 1) * variance_cell[j] + (n2 * n2) / (n2 + n2)
                                   * (mean_cell[i] - mean_cell[j]) ** 2) / (n2 + n2 - 1))

            if combined_variance1 > combined_variance2:
                larger_variance = combined_variance1
                smaller_variance = combined_variance2
                df_larger = n1 + n1 - 1
                df_smaller = n2 + n2 - 1
            else:
                larger_variance = combined_variance2
                smaller_variance = combined_variance1
                df_larger = n2 + n2 - 1
                df_smaller = n1 + n1 - 1

            if smaller_variance > 0:
                f_statistic = larger_variance / smaller_variance
                f_critical = f.ppf(1 - p_value, df_larger, df_smaller)
                csn[i, j] = f_statistic > f_critical

    return csr_matrix(csn)


def difference_csr(A: csr_matrix, B: csr_matrix) -> csr_matrix:
    # Convert non-zero elements to 1
    A.data[A.data > 0] = 1
    B.data[B.data > 0] = 1
    # Subtract B from A
    C = A - B
    C.data[C.data < 0] = 0
    C.eliminate_zeros()
    return C


def csnet(
        adata,
        mean_all,
        variance_all,
        cell_num_all,
        cell_name,
        total_gene_dict,
        rt_matrix,
        tf_matrix,
        alpha: float = 0.01,
):
    data = adata.X.T.toarray()
    gene_list = adata.var_names
    n1, n2 = data.shape  # n1是基因数，n2是细胞数
    n3 = len(total_gene_dict)  # n3是总基因数

    csn1 = pearson_net(data, alpha)
    mean_cell = np.mean(data, axis=1)  # Calculate the mean of selected rows
    variance_cell = np.var(data, axis=1)  # Calculate the variance of selected rows

    shape = (n3, n3)
    rp_tf_rows = []
    rp_tf_cols = []
    rp_tf_data = []
    tf_tg_rows = []
    tf_tg_cols = []
    tf_tg_data = []

    if rt_matrix is not None and rt_matrix.nnz > 0:
        dif1 = difference_csr(rt_matrix, csn1)
        pearson_rt = csn1.multiply(rt_matrix)
        f_rt = variance_direct(dif1, gene_list, mean_all, variance_all, cell_num_all, mean_cell, variance_cell, n2, alpha)
        combined = pearson_rt + f_rt

        rows, cols = combined.nonzero()
        for i, j in zip(rows, cols):
            gene1, gene2 = gene_list[i], gene_list[j]
            gene1, gene2 = total_gene_dict[gene1], total_gene_dict[gene2]
            rp_tf_rows.append(gene1)
            rp_tf_cols.append(gene2)
            rp_tf_data.append(1)

    if tf_matrix is not None and tf_matrix.nnz > 0:
        dif2 = difference_csr(tf_matrix, csn1)
        pearson_tf1 = csn1.multiply(tf_matrix)
        f_tf1 = variance_direct(dif2, gene_list, mean_all, variance_all, cell_num_all, mean_cell, variance_cell, n2, alpha)
        combined = pearson_tf1 + f_tf1

        rows, cols = combined.nonzero()
        for i, j in zip(rows, cols):
            gene1, gene2 = gene_list[i], gene_list[j]
            gene1, gene2 = total_gene_dict[gene1], total_gene_dict[gene2]
            tf_tg_rows.append(gene1)
            tf_tg_cols.append(gene2)
            tf_tg_data.append(1)

    rp_tf_adj = coo_matrix((rp_tf_data, (rp_tf_rows, rp_tf_cols)), shape=shape)
    tf_tg_adj = coo_matrix((tf_tg_data, (tf_tg_rows, tf_tg_cols)), shape=shape)
    print('Cell type', cell_name, 'is completed')

    return rp_tf_adj, tf_tg_adj


    '''
        if len(rp_tf) > 0:
        os.makedirs('Receptor_TF_network', exist_ok=True)
        with open(f"Receptor_TF_network/{cell_name}.txt", 'w') as f:
            f.write(f"Receptor\tTF\tType\n")
            for key, value in rp_tf.items():
                gene1 = gene_dict[key]
                for tf in value:
                    gene2 = gene_dict[tf]
                    typ = interaction_type[gene1][gene2]
                    f.write(f"{key}\t{tf}\t{typ}\n")

    if len(tf_tg) > 0:
        os.makedirs('TF_TG_network', exist_ok=True)
        with open(f"TF_TG_network/{cell_name}.txt", 'w') as f:
            f.write(f"TF\tTG\tType\n")
            for key, value in tf_tg.items():
                for tg in value:
                    f.write(f"{key}\t{tg}\tTF-target\n")
    '''


