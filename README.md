# CellSigN

CellSigN reconstructs candidate **ligand–receptor–mediator–transcription
factor–target (L–R–M–TF–TG)** signaling paths from single-cell RNA-seq data.
It combines receiver-cell statistical evidence with curated ligand–receptor
and intracellular references. By default, it analyzes both directions between
the malignant cell type and every other eligible cell type; specific ordered
sender–receiver pairs can also be requested.

**Inputs:** an `.h5ad` file or a gene-by-cell expression table with cell labels,
plus the two bundled reference tables. **Outputs:** candidate pathway tables,
receiver-edge evidence, DEG lists and a run manifest. The source module is
run with `python -m main`; installation also provides the `cellsign` command.

## Quick start

Use Python 3.9 or newer. Extract the repository archive and run these commands
from its top-level directory. In an existing compatible environment, skip the
virtual-environment creation and activation commands.

Windows CMD:

```cmd
python -m venv .venv
.venv\Scripts\activate
python -m pip install -e .
```

Linux/macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The following **small synthetic smoke test** uses files in `tests/data/` and
works as one line in Windows CMD, PowerShell or Bash:

```text
python -m main --stage all --expression tests/data/expression.tsv --cell-path tests/data/cell_types.tsv --intracellular-prior tests/data/intracellular_prior.tsv --ligand-receptor-prior tests/data/ligand_receptor.tsv --pair Sender:Receiver --stability-resamples 0 --path-permutations 5 --output-dir Results_smoke
```

A successful run ends with `CellSigN stage 'all' completed.` and writes
`Results_smoke_expression/Sender_to_Receiver_pathway.txt`. This synthetic
fixture checks installation and execution; it is not a biological case study.
To run the automated test suite, install `.[test]` and execute `pytest`.

## Run the included Chung breast cancer example

`Example data/` contains a processed **515-cell × 5,000-gene** expression
matrix and matching cell labels derived from the primary breast cancer dataset
of Chung et al. (2017), available from GEO as
[GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688).
The table exports the `X` matrix of the processed h5ad used for analysis;
it is **not** the original GEO raw matrix. See [data provenance](DATA_PROVENANCE.md)
for details.

Run the supplied text export with the two bundled reference tables:

```text
python -m main --stage all --expression "Example data/Data_Chung2017_Breast_all_expression.txt" --cell-path "Example data/Data_Chung2017_Breast_all_cell_labels.txt" --intracellular-prior reference_data/human/intracellular_network.txt.gz --ligand-receptor-prior reference_data/human/ligand_receptor.txt --malignant-label Malignant --no-normalize --no-log1p --de-threshold 0.05 --k-paths 1 --path-permutations 500 --seed 42 --output-dir Results
```

The results are written under `Results_Data_Chung2017_Breast_all_expression/`.
This real-data run is more computationally demanding than the smoke test.
`Example Result/` contains representative B-cell-to-malignant and
malignant-to-B-cell pathway tables from a run on the processed h5ad. The text
export may yield small numerical differences and is not guaranteed to
reproduce every row exactly. These tables are candidate associations, not
experimentally validated signaling paths.

## Use your own data

- For `.h5ad`, cells are rows and genes are columns. Cell labels must be in
  `adata.obs`; the default column is `celltype`, configurable with
  `--cell-type-column`.
- For a delimited text matrix, genes are rows, cells are columns and the first
  column contains gene identifiers. Supply `--cell-path` with a two-column
  cell-ID/cell-type file; its first row may be a header.
- Normalization and `log1p` are enabled by default. Use `--no-normalize` and
  `--no-log1p` when the input values are already processed and should be used
  unchanged.

The intracellular reference is bundled as
`reference_data/human/intracellular_network.txt.gz`. CellSigN reads it
directly; decompression is unnecessary. The ligand–receptor reader uses the
first two columns of `reference_data/human/ligand_receptor.txt` and does not
apply a score threshold because this supplied table is already curated.
Reference origins, limitations and checksums are described in
[data provenance](DATA_PROVENANCE.md) and
[source provenance](SOURCE_PROVENANCE.md).

## How CellSigN works

1. Filter cells and genes, optionally normalize, and select up to 5,000 HVGs
   across all cells.
2. For each eligible cell type, run one-versus-rest Wilcoxon differential
   expression once. Keep all HVGs with positive logFC and raw p-value below
   `--de-threshold`; there is no post-DE Top-N cutoff.
3. Within each receiver cell type, test every unordered pair of receiver DEGs
   using Pearson correlation and a nested linear-versus-cubic model comparison.
   Apply BH correction to the edge-test families and assess effect size and
   bootstrap stability.
4. Retain evidence-supported association edges and assign edge costs. Intersect
   these edges with the intracellular signaling reference for path search;
   reference directions and receptor/TF roles constrain the resulting graph.
   A prior edge without statistical support does not enter the path graph.
5. Find up to `--k-paths` lowest-cost simple receptor-to-TF paths, assess
   path and receptor costs by edge-cost permutation, and attach sender ligands
   and receiver TF targets. Ordered directions are evaluated independently.

The [command-line guide](README_COMMAND_LINE_GUIDE.md)
shows one-stage and separate `evidence`/`paths` runs. A `paths` run reuses
saved evidence rather than recalculating receiver statistics.

## Read the results

CellSigN appends the expression filename to `--output-dir`, creating a folder
such as `Results_Data_Chung2017_Breast_all_expression/`. Its key contents are:

| File or directory | Contents |
|---|---|
| `<Sender>_to_<Receiver>_pathway.txt` | Candidate L–R–M–TF–TG pathways for one ordered cell direction |
| `gene_lists/` | Selected one-versus-rest DEGs by cell type |
| `evidence/` | Full DE and receiver-edge evidence, path edges and R–M–TF subnetworks |
| `other/differential_gene_counts.txt` | One DEG count per eligible cell type |
| `other/test_counts.txt` | Per-direction counts of tested pairs and permutations |
| `other/run_manifest.json` | Input hashes, parameters, software versions and output counts |

The first five pathway columns are `Ligand`, `Receptor`, `Mediator`, `TF` and
`Target`. **Each row represents one intracellular R–M–TF path.** Commas in
`Mediator` list sequential nodes in path order; a direct R–TF path has no
mediator. Commas in `Ligand` or `Target` list alternative partners attached to
that same intracellular path.

`Path_P`/`Path_Q` are R–TF pair-level permutation results;
`Receptor_P`/`Receptor_Q` summarize the best reachable TF for a receptor.
`Receptor_Significant` compares `Receptor_Q` with `--path-alpha`. These
statistics are reported, but nonsignificant pathways are not automatically
deleted. See [METHODS.md](METHODS.md) before interpreting their test families.

## Main options and reproducibility

| Option | Default | Effect |
|---|---:|---|
| `--stage` | `all` | Run both stages, only `evidence`, or only `paths` |
| `--malignant-label` | `Malignant` | Label used for the default malignant↔other scope |
| `--pair` | none | Analyze one ordered `SENDER:RECEIVER` pair; repeat to select more |
| `--hvg-top-genes` | 5000 | Maximum number of globally selected HVGs |
| `--de-threshold` | 0.05 | Raw one-versus-rest Wilcoxon p-value cutoff |
| `--alpha` | 0.05 | BH-adjusted Pearson and nonlinear edge-test cutoff |
| `--stability-resamples` | 20 | Receiver-cell bootstrap iterations |
| `--k-paths` | 5 | Lowest-cost paths retained per R–TF pair |
| `--path-permutations` | 500 | Edge-cost permutations; `0` disables them |
| `--path-alpha` | 0.05 | Receptor-q threshold for the significance flag |
| `--seed` | 0 | Random seed |

Use `python -m main --help` for all options, including the cost weights.
Cell types with fewer than five cells are excluded. The run manifest records
the exact inputs, hashes, settings and software versions needed to assess
reproducibility. CellSigN reports **association-supported candidate paths**;
differential expression, statistical association and permutation tests do not
by themselves establish causal signaling. Independent biological validation
is required for causal claims.

The source code is distributed under [GPL-3.0](LICENSE). The bundled reference
tables may also be subject to their upstream sources' terms; see
[data provenance](DATA_PROVENANCE.md).
