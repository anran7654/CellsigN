"""Tests for CellSigN statistical utilities."""

import numpy as np

from main.statistics import (
    benjamini_hochberg,
    cubic_nonlinear_receiver_pair_tests,
    fisher_correlation_difference,
    pearson_matrix,
    two_group_anova,
)


def test_benjamini_hochberg_known_values():
    p = np.array([0.01, 0.04, 0.03, 0.002])
    q = benjamini_hochberg(p)
    np.testing.assert_allclose(q, [0.02, 0.04, 0.04, 0.008])


def test_benjamini_hochberg_preserves_nan():
    q = benjamini_hochberg(np.array([0.01, np.nan, 0.5]))
    assert np.isnan(q[1])
    np.testing.assert_allclose(q[[0, 2]], [0.02, 0.5])


def test_two_group_anova_detects_shift_and_reports_effect_size():
    receiver = np.column_stack((np.arange(10.0) + 10, np.arange(10.0)))
    background = np.column_stack((np.arange(10.0), np.arange(10.0)))
    statistic, p_value, eta_squared = two_group_anova(receiver, background)
    assert statistic[0] > 10 and p_value[0] < 0.01 and eta_squared[0] > 0.4
    assert statistic[1] == 0 and p_value[1] == 1 and eta_squared[1] == 0


def test_pearson_and_fisher_difference_use_disjoint_group_sizes():
    x = np.arange(10.0)
    receiver = np.column_stack((x, x))
    background = np.column_stack((x, x[::-1]))
    receiver_r, _ = pearson_matrix(receiver)
    background_r, _ = pearson_matrix(background)
    z, p = fisher_correlation_difference(
        receiver_r[0, 1], len(receiver), background_r[0, 1], len(background)
    )
    assert z > 0 and p < 0.001


def test_receiver_only_nested_anova_detects_cubic_increment():
    x = np.linspace(-2, 2, 30)
    expression = np.column_stack((x, x**3))
    result = cubic_nonlinear_receiver_pair_tests(
        expression, np.array([0]), np.array([1])
    )
    _, p_uv, delta_uv, _, _, _ = result
    assert p_uv[0] < 0.001
    assert delta_uv[0] > 0.5
