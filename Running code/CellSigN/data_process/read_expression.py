import csv

import numpy as np


def read_gene_expression(file_path):
    sample_names = []  # One-dimensional vector for sample names
    gene_names = []  # One-dimensional vector for gene names
    expression_matrix = []  # Two-dimensional array for expression values

    with open(file_path, 'r', encoding='utf-8') as file:
        line1 = next(file)
        sample_names = line1.strip().replace(',', '\t').split('\t')
        for line in file:
            parts = line.strip().replace(',', '\t').split('\t')
            gene_names.append(parts[0])
            expression_matrix.append([convert2float(value) for value in parts[1:]])
    data = np.array(expression_matrix).T

    return sample_names, gene_names, data


def convert2float(value):
    try:
        # Attempt to convert the value to a float
        float_value = float(value)
        return float_value
    except ValueError:
        # If conversion fails, output 0
        return 0
