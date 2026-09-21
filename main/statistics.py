"""Statistical utilities used by CellSigN.

The functions in this module operate on disjoint receiver and background
cells. They return raw statistics and p-values; the caller applies correction
after every hypothesis in a cell type has been evaluated.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import f, norm, t


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values, preserving NaNs."""
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return adjusted
    p = np.clip(values[finite], 0.0, 1.0)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    n = ranked.size
    q_ranked = ranked * n / np.arange(1, n + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0.0, 1.0)
    q = np.empty_like(q_ranked)
    q[order] = q_ranked
    adjusted[finite] = q
    return adjusted


def q_to_score(q_values: np.ndarray, epsilon: float = 1e-300) -> np.ndarray:
    """Map adjusted p-values to a bounded evidence score in [0, 1)."""
    q = np.asarray(q_values, dtype=float)
    score = np.zeros(q.shape, dtype=float)
    finite = np.isfinite(q)
    if finite.any():
        clipped = np.clip(q[finite], epsilon, 1.0)
        strength = -np.log10(clipped)
        score[finite] = strength / (1.0 + strength)
    return score


def p_to_score(p_values: np.ndarray, epsilon: float = 1e-300) -> np.ndarray:
    """Map raw p-values to a bounded evidence score in [0, 1)."""
    p = np.asarray(p_values, dtype=float)
    score = np.zeros(p.shape, dtype=float)
    finite = np.isfinite(p)
    if finite.any():
        clipped = np.clip(p[finite], epsilon, 1.0)
        strength = -np.log10(clipped)
        score[finite] = strength / (1.0 + strength)
    return score


def two_group_anova(
    receiver: np.ndarray,
    background: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return one-way two-group F, p, and eta-squared for every gene."""
    receiver = np.asarray(receiver, dtype=float)
    background = np.asarray(background, dtype=float)
    if receiver.ndim != 2 or background.ndim != 2:
        raise ValueError("ANOVA inputs must be two-dimensional matrices.")
    if receiver.shape[1] != background.shape[1]:
        raise ValueError("Receiver and background matrices must have the same genes.")
    n1, n0 = receiver.shape[0], background.shape[0]
    if n1 < 2 or n0 < 2:
        raise ValueError("ANOVA requires at least two receiver and two background cells.")
    mean1 = receiver.mean(axis=0)
    mean0 = background.mean(axis=0)
    grand = (n1 * mean1 + n0 * mean0) / (n1 + n0)
    ss_between = n1 * (mean1 - grand) ** 2 + n0 * (mean0 - grand) ** 2
    ss_within = ((receiver - mean1) ** 2).sum(axis=0) + ((background - mean0) ** 2).sum(axis=0)
    df_within = n1 + n0 - 2
    ms_within = ss_within / df_within
    f_stat = np.zeros(receiver.shape[1], dtype=float)
    positive_error = ms_within > 0
    f_stat[positive_error] = ss_between[positive_error] / ms_within[positive_error]
    zero_error_difference = (~positive_error) & (ss_between > 0)
    f_stat[zero_error_difference] = np.inf
    p_value = f.sf(f_stat, 1, df_within)
    total = ss_between + ss_within
    eta_squared = np.divide(
        ss_between, total, out=np.zeros_like(ss_between), where=total > 0
    )
    return f_stat, p_value, np.clip(eta_squared, 0.0, 1.0)


def pearson_matrix(expression: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return pairwise Pearson coefficients and two-sided p-values."""
    expression = np.asarray(expression, dtype=float)
    if expression.ndim != 2:
        raise ValueError("Expression must be cells by genes.")
    n = expression.shape[0]
    if n < 3:
        raise ValueError("Pearson tests require at least three cells.")
    centered = expression - expression.mean(axis=0, keepdims=True)
    norms = np.sqrt((centered**2).sum(axis=0))
    denominator = np.outer(norms, norms)
    r = np.divide(
        centered.T @ centered,
        denominator,
        out=np.zeros((expression.shape[1], expression.shape[1]), dtype=float),
        where=denominator > 0,
    )
    r = np.clip(r, -1.0, 1.0)
    df = n - 2
    residual = np.maximum(1.0 - r**2, np.finfo(float).tiny)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        statistic = r * np.sqrt(df / residual)
    p_value = 2.0 * t.sf(np.abs(statistic), df)
    constant = norms == 0
    if constant.any():
        r[constant, :] = 0.0
        r[:, constant] = 0.0
        p_value[constant, :] = 1.0
        p_value[:, constant] = 1.0
    np.fill_diagonal(r, np.nan)
    np.fill_diagonal(p_value, np.nan)
    return r, p_value


def pearson_pairs(
    expression: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Pearson coefficients and p-values for selected gene pairs only."""
    expression = np.asarray(expression, dtype=float)
    source_index = np.asarray(source_index, dtype=int)
    target_index = np.asarray(target_index, dtype=int)
    if expression.ndim != 2 or source_index.shape != target_index.shape:
        raise ValueError("Pearson pair inputs have incompatible shapes.")
    n = expression.shape[0]
    if n < 3:
        raise ValueError("Pearson tests require at least three cells.")
    x = expression[:, source_index]
    y = expression[:, target_index]
    x = x - x.mean(axis=0, keepdims=True)
    y = y - y.mean(axis=0, keepdims=True)
    denominator = np.sqrt((x * x).sum(axis=0) * (y * y).sum(axis=0))
    r = np.divide(
        (x * y).sum(axis=0), denominator,
        out=np.zeros(source_index.size, dtype=float), where=denominator > 0,
    )
    r = np.clip(r, -1.0, 1.0)
    residual = np.maximum(1.0 - r**2, np.finfo(float).tiny)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        statistic = r * np.sqrt((n - 2) / residual)
    p_value = 2.0 * t.sf(np.abs(statistic), n - 2)
    p_value[denominator == 0] = 1.0
    return r, p_value


def fisher_correlation_difference(
    receiver_r: np.ndarray,
    receiver_n: int,
    background_r: np.ndarray,
    background_n: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compare independent receiver and background Pearson correlations."""
    receiver_r, background_r = np.broadcast_arrays(
        np.asarray(receiver_r, dtype=float), np.asarray(background_r, dtype=float)
    )
    z_stat = np.full(receiver_r.shape, np.nan, dtype=float)
    p_value = np.full(receiver_r.shape, np.nan, dtype=float)
    if receiver_n <= 3 or background_n <= 3:
        return z_stat, p_value
    valid = np.isfinite(receiver_r) & np.isfinite(background_r)
    if not valid.any():
        return z_stat, p_value
    # Sample correlations numerically indistinguishable from +/-1 otherwise
    # produce unstable differences after arctanh.
    eps = 1e-7
    zr = np.arctanh(np.clip(receiver_r[valid], -1.0 + eps, 1.0 - eps))
    zb = np.arctanh(np.clip(background_r[valid], -1.0 + eps, 1.0 - eps))
    standard_error = np.sqrt(1.0 / (receiver_n - 3) + 1.0 / (background_n - 3))
    z_stat[valid] = (zr - zb) / standard_error
    p_value[valid] = 2.0 * norm.sf(np.abs(z_stat[valid]))
    return z_stat, p_value


def _rss(design: np.ndarray, outcomes: np.ndarray) -> tuple[np.ndarray, int, int]:
    coefficients, _, rank, _ = np.linalg.lstsq(design, outcomes, rcond=None)
    residual = outcomes - design @ coefficients
    return (residual**2).sum(axis=0), int(rank), int(design.shape[0] - rank)


def _standardized_predictor(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    scale = np.sqrt(np.mean(centered**2))
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(centered)
    return centered / scale


def cubic_nonlinear_tests(
    expression: np.ndarray,
    receiver_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit directional cubic nested models for every ordered gene pair.

    The reduced model contains group-specific linear effects and common
    quadratic/cubic terms. The full model adds receiver-specific quadratic
    and cubic interactions. Rows index predictors and columns index outcomes.
    """
    expression = np.asarray(expression, dtype=float)
    receiver_mask = np.asarray(receiver_mask, dtype=bool)
    if expression.ndim != 2 or receiver_mask.shape != (expression.shape[0],):
        raise ValueError("Nonlinear-test inputs have incompatible shapes.")
    if receiver_mask.sum() < 4 or (~receiver_mask).sum() < 2:
        size = expression.shape[1]
        empty = np.full((size, size), np.nan, dtype=float)
        return empty.copy(), empty.copy(), empty.copy()
    n, genes = expression.shape
    group = receiver_mask.astype(float)
    outcomes = expression
    f_matrix = np.full((genes, genes), np.nan, dtype=float)
    p_matrix = np.full((genes, genes), np.nan, dtype=float)
    delta_matrix = np.full((genes, genes), np.nan, dtype=float)
    receiver_outcomes = outcomes[receiver_mask]
    for predictor_index in range(genes):
        z = _standardized_predictor(expression[:, predictor_index])
        if np.all(z == 0):
            continue
        z2 = z**2
        z3 = z**3
        reduced = np.column_stack((np.ones(n), group, z, group * z, z2, z3))
        full = np.column_stack((reduced, group * z2, group * z3))
        rss_reduced, rank_reduced, _ = _rss(reduced, outcomes)
        rss_full, rank_full, df_full = _rss(full, outcomes)
        numerator_df = rank_full - rank_reduced
        if numerator_df > 0 and df_full > 0:
            improvement = np.maximum(rss_reduced - rss_full, 0.0)
            denominator = rss_full / df_full
            f_stat = np.divide(
                improvement / numerator_df,
                denominator,
                out=np.zeros_like(improvement),
                where=denominator > 0,
            )
            f_stat[(denominator == 0) & (improvement > 0)] = np.inf
            f_matrix[predictor_index] = f_stat
            p_matrix[predictor_index] = f.sf(f_stat, numerator_df, df_full)
        zr = _standardized_predictor(expression[receiver_mask, predictor_index])
        if not np.all(zr == 0):
            linear = np.column_stack((np.ones(zr.size), zr))
            cubic = np.column_stack((linear, zr**2, zr**3))
            rss_linear, _, _ = _rss(linear, receiver_outcomes)
            rss_cubic, _, _ = _rss(cubic, receiver_outcomes)
            delta = np.divide(
                np.maximum(rss_linear - rss_cubic, 0.0),
                rss_linear,
                out=np.zeros_like(rss_linear),
                where=rss_linear > 0,
            )
            delta_matrix[predictor_index] = np.clip(delta, 0.0, 1.0)
    np.fill_diagonal(f_matrix, np.nan)
    np.fill_diagonal(p_matrix, np.nan)
    np.fill_diagonal(delta_matrix, np.nan)
    return f_matrix, p_matrix, delta_matrix


def cubic_increment_effect(x: np.ndarray, y: np.ndarray) -> float:
    """Return target-only cubic-versus-linear incremental R-squared."""
    x = _standardized_predictor(np.asarray(x, dtype=float))
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    if x.size < 4 or np.all(x == 0):
        return 0.0
    linear = np.column_stack((np.ones(x.size), x))
    cubic = np.column_stack((linear, x**2, x**3))
    rss_linear, _, _ = _rss(linear, y)
    rss_cubic, _, _ = _rss(cubic, y)
    if rss_linear[0] <= 0:
        return 0.0
    return float(np.clip((rss_linear[0] - rss_cubic[0]) / rss_linear[0], 0.0, 1.0))


def cubic_nonlinear_receiver_pair_tests(
    receiver_expression: np.ndarray,
    source_index: np.ndarray,
    target_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Test cubic-over-linear increments within receiver cells for gene pairs.

    Each unordered pair is fitted in both predictive directions.  The reduced
    model contains an intercept and linear term; the full model adds quadratic
    and cubic terms.  The returned delta is the incremental R-squared.
    """
    expression = np.asarray(receiver_expression, dtype=float)
    source_index = np.asarray(source_index, dtype=int)
    target_index = np.asarray(target_index, dtype=int)
    if expression.ndim != 2:
        raise ValueError("Receiver expression must be a cells-by-genes matrix.")
    if source_index.shape != target_index.shape:
        raise ValueError("Nonlinear pair indices have incompatible shapes.")
    size = source_index.size
    f_uv = np.full(size, np.nan, dtype=float)
    p_uv = np.full(size, np.nan, dtype=float)
    d_uv = np.full(size, np.nan, dtype=float)
    f_vu = np.full(size, np.nan, dtype=float)
    p_vu = np.full(size, np.nan, dtype=float)
    d_vu = np.full(size, np.nan, dtype=float)
    if expression.shape[0] < 5 or size == 0:
        return f_uv, p_uv, d_uv, f_vu, p_vu, d_vu

    def fit_direction(predictor: int, outcome: int) -> tuple[float, float, float]:
        n = expression.shape[0]
        z = _standardized_predictor(expression[:, predictor])
        if np.all(z == 0):
            return np.nan, np.nan, 0.0
        y = expression[:, outcome:outcome + 1]
        reduced = np.column_stack((np.ones(n), z))
        full = np.column_stack((reduced, z**2, z**3))
        rss_reduced, rank_reduced, _ = _rss(reduced, y)
        rss_full, rank_full, df_full = _rss(full, y)
        numerator_df = rank_full - rank_reduced
        f_stat = np.nan
        p_value = np.nan
        if numerator_df > 0 and df_full > 0:
            improvement = max(float(rss_reduced[0] - rss_full[0]), 0.0)
            denominator = float(rss_full[0] / df_full)
            if denominator > 0:
                f_stat = (improvement / numerator_df) / denominator
            elif improvement > 0:
                f_stat = np.inf
            else:
                f_stat = 0.0
            p_value = float(f.sf(f_stat, numerator_df, df_full))
        delta = 0.0
        if rss_reduced[0] > 0:
            delta = float(np.clip(
                (rss_reduced[0] - rss_full[0]) / rss_reduced[0], 0.0, 1.0
            ))
        return float(f_stat), float(p_value), float(delta)

    for index, (u, v) in enumerate(zip(source_index, target_index)):
        f_uv[index], p_uv[index], d_uv[index] = fit_direction(int(u), int(v))
        f_vu[index], p_vu[index], d_vu[index] = fit_direction(int(v), int(u))
    return f_uv, p_uv, d_uv, f_vu, p_vu, d_vu
