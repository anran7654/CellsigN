import numpy as np
from scipy.stats import t, f, rankdata
from scipy.sparse import csr_matrix


def pearson_correlation(data):
    data_centered = data - np.mean(data, axis=1, keepdims=True)
    std_devs = np.std(data_centered, axis=1, keepdims=True)
    zero_std_mask = (std_devs == 0)  # Identify where the standard deviations are zero
    std_devs[zero_std_mask] = 1  # Replace zero standard deviations with ones to avoid division by zero
    data_normalized = data_centered / std_devs
    matrix = np.dot(data_normalized, data_normalized.T) / (data.shape[1] - 1)  # Compute the correlation matrix
    np.fill_diagonal(matrix, 0)  # setting all the elements on the main diagonal to zero
    return matrix


def spearman_correlation(data):
    ranked_data = np.apply_along_axis(rankdata, 1, data)  # Rank the data
    return pearson_correlation(ranked_data)  # Use the fast Pearson correlation on ranked data


def var_net(variance_all, n1, variance_cell, n2, alpha=0.01):
    n3 = len(variance_all)

    # Find the larger and smaller elements at each position
    larger = np.maximum(variance_all, variance_cell)
    smaller = np.minimum(variance_all, variance_cell)

    # Calculate the ratio
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = np.divide(larger, smaller)
        ratios = np.where(np.isfinite(ratios), ratios, 0)

    # Determine which array element is larger at each position
    array1_is_larger = variance_all > variance_cell
    df_larger = np.where(array1_is_larger, n1, n2)
    df_smaller = np.where(array1_is_larger, n2, n1)

    f_critical = f.ppf(1 - alpha, df_larger, df_smaller)

    sig = ratios > f_critical
    csn = np.where(sig, 1, 0)
    return csn


def statistic(matrix, n_samples, p_value=0.01):
    epsilon = 1e-10  # Clipping the correlation values to avoid invalid values in sqrt
    matrix = np.clip(matrix, -1 + epsilon, 1 - epsilon)
    df = n_samples - 2  # Degree of freedom
    t_stat = matrix * np.sqrt(df / (1 - matrix ** 2))  # t-statistics
    # pval_matrix = 2 * (1 - t.cdf(np.abs(t_stat), df))
    alpha = p_value / 2
    pt = t.ppf(1 - alpha, df)
    csn = csr_matrix(np.abs(t_stat) > pt)
    return csn
