import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from main.io import (
    one_vs_rest_cell_type_de_tables,
    select_cell_type_gene_tables,
)


def test_top_gene_ranking_skips_cell_types_too_small_for_network_statistics():
    expression = np.log1p(np.array([
        [9, 1, 1], [8, 1, 0], [7, 2, 1], [9, 1, 2], [8, 2, 1],
        [1, 9, 1], [1, 8, 0], [2, 7, 1], [1, 9, 2], [2, 8, 1],
        [1, 1, 9],
    ], dtype=float))
    adata = AnnData(
        X=expression,
        obs=pd.DataFrame({
            "celltype": pd.Categorical(["A"] * 5 + ["B"] * 5 + ["singleton"])
        }),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )

    with pytest.warns(RuntimeWarning, match="singleton"):
        tables = select_cell_type_gene_tables(adata, "celltype", 2)

    assert list(tables) == ["A", "B"]


def test_top_gene_ranking_rejects_data_without_an_eligible_cell_type():
    adata = AnnData(
        X=np.ones((4, 2), dtype=float),
        obs=pd.DataFrame({"celltype": pd.Categorical(["A", "A", "B", "B"])}),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    with pytest.raises(ValueError, match="No cell type has at least five"):
        select_cell_type_gene_tables(adata, "celltype", 2)


def test_one_vs_rest_de_is_computed_once_per_requested_cell_type():
    expression = np.log1p(np.array(
        [[12, 1, 1], [11, 1, 1], [10, 2, 1], [13, 1, 2], [12, 2, 1]]
        + [[1, 12, 1], [1, 11, 1], [2, 10, 1], [1, 13, 2], [2, 12, 1]]
        + [[1, 1, 12], [1, 1, 11], [1, 2, 10], [2, 1, 13], [1, 2, 12]],
        dtype=float,
    ))
    adata = AnnData(
        X=expression,
        obs=pd.DataFrame({
            "celltype": pd.Categorical(["A"] * 5 + ["B"] * 5 + ["C"] * 5)
        }),
        var=pd.DataFrame(index=["GA", "GB", "GC"]),
    )
    raw = one_vs_rest_cell_type_de_tables(adata, "celltype", ["A", "C"])
    assert set(raw) == {"A", "C"}
    assert raw["A"]["Reference"].eq("All_other_cells").all()
    assert "P_Value" in raw["A"].columns
    assert "Q_Value" not in raw["A"].columns
