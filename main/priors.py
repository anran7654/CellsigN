"""Readers for intracellular and ligand-receptor prior networks."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


TF_TG_TYPE = "controls-expression-of"
REVERSIBLE_TYPES = {"interacts-with", "in-complex-with"}


def load_intracellular_prior(
    path: str | Path,
    genes: set[str],
    chunksize: int = 250_000,
) -> pd.DataFrame:
    """Load prior edges whose endpoints occur in the selected gene universe."""
    columns = ["source", "target", "prior_type", "layer"]
    kept: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        chunksize=chunksize,
        dtype=str,
    ):
        required = {"Gene1", "Gene2", "Type"}
        if not required.issubset(chunk.columns):
            raise ValueError(f"Intracellular prior must contain columns {sorted(required)}.")
        chunk = chunk.loc[
            chunk["Gene1"].isin(genes) & chunk["Gene2"].isin(genes),
            ["Gene1", "Gene2", "Type"],
        ].rename(columns={"Gene1": "source", "Gene2": "target", "Type": "prior_type"})
        if chunk.empty:
            continue
        chunk["layer"] = chunk["prior_type"].eq(TF_TG_TYPE).map(
            {True: "tf_target", False: "signaling"}
        )
        kept.append(chunk[columns])
        reverse = chunk[chunk["prior_type"].isin(REVERSIBLE_TYPES)].copy()
        if not reverse.empty:
            reverse[["source", "target"]] = reverse[["target", "source"]].to_numpy()
            kept.append(reverse[columns])
    if not kept:
        return pd.DataFrame(columns=columns)
    data = pd.concat(kept, ignore_index=True).drop_duplicates()
    data = data[data["source"] != data["target"]]
    data = (
        data.groupby(["source", "target", "layer"], as_index=False, sort=False)["prior_type"]
        .agg(lambda values: ";".join(sorted(set(values))))
    )
    return data[columns]


def load_ligand_receptor_prior(
    path: str | Path,
) -> pd.DataFrame:
    """Load every ligand-receptor record using only the first two columns.

    The bundled table is already curated. Any additional columns are ignored,
    and no score threshold is applied. Expression-universe filtering happens
    later for each ordered sender-receiver analysis.
    """
    table = pd.read_csv(
        path, sep="\t", compression="infer", usecols=[0, 1], dtype=str
    )
    if table.shape[1] != 2:
        raise ValueError("Ligand-receptor prior must contain at least two columns.")
    table.columns = ["Ligand", "Receptor"]
    table = table.dropna(subset=["Ligand", "Receptor"]).copy()
    table["Ligand"] = table["Ligand"].str.strip()
    table["Receptor"] = table["Receptor"].str.strip()
    table = table[(table["Ligand"] != "") & (table["Receptor"] != "")]
    return table.reset_index(drop=True)
