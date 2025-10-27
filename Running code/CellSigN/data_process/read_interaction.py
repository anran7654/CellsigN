import os
import numpy as np
from scipy.sparse import csr_matrix


def pathway(directory, gene_dict, col1_index: int = 0, col2_index: int = 1, col3_index: int = 2):
    n1 = len(gene_dict)
    rp_tf = np.zeros((n1, n1), dtype=int)
    tf_tg = np.zeros((n1, n1), dtype=int)
    interaction_type = {}
    if os.path.isfile(directory):  # Check if it's a file (not a directory)
        with open(directory, 'r', encoding='utf-8') as file:
            next(file)
            for line in file:
                elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
                gene1 = elements[col1_index]
                gene2 = elements[col2_index]
                c3 = elements[col3_index]
                if gene1 in gene_dict and gene2 in gene_dict:
                    idx1, idx2 = gene_dict[gene1], gene_dict[gene2]
                    if c3 == 'controls-expression-of':
                        tf_tg[idx1, idx2] = 1
                    else:
                        rp_tf[idx1, idx2] = 1

                    pair = gene1+'_'+gene2
                    if pair not in interaction_type:
                        interaction_type[pair] = set()
                        interaction_type[pair].add(c3)
                    else:
                        interaction_type[pair].add(c3)

                    if c3 == 'interacts-with' or c3 == 'in-complex-with':
                        rp_tf[idx2, idx1] = 1

                        pair = gene2 + '_' + gene1
                        if pair not in interaction_type:
                            interaction_type[pair] = set()
                            interaction_type[pair].add(c3)
                        else:
                            interaction_type[pair].add(c3)

        np.fill_diagonal(tf_tg, 0)
        np.fill_diagonal(rp_tf, 0)

        return csr_matrix(rp_tf), csr_matrix(tf_tg), interaction_type
    else:
        print("Error: The reference network file does not exist.")
        return csr_matrix(rp_tf), csr_matrix(tf_tg), interaction_type


def ligand_receptor(directory, gene_dict, col1_index: int = 0, col2_index: int = 1, col3_index: int = 2, min_score: int = 1):
    lg_rp = {}
    if os.path.isfile(directory):  # Check if it's a file (not a directory)
        with open(directory, 'r', encoding='utf-8') as file:
            next(file)
            for line in file:
                elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
                gene1 = elements[col1_index]
                gene2 = elements[col2_index]
                score = elements[col3_index]
                if int(score) >= min_score and gene1 in gene_dict and gene2 in gene_dict:
                    idx1, idx2 = gene_dict[gene1], gene_dict[gene2]
                    if idx1 not in lg_rp:
                        lg_rp[idx1] = set()  # Create a new set if it doesn't exist
                    lg_rp[idx1].add(idx2)
        return lg_rp
    else:
        print("Error: The ligand-receptor file does not exist.")
        return lg_rp
