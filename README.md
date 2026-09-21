# CellSigN

The installable distribution and optional console command are named
`cellsign`; the Python module entry point is `python -m main`.

CellSigN infers candidate bidirectional `Ligand-Receptor-Mediator-TF-Target`
paths from single-cell RNA-seq data. The current workflow tests each eligible
cell type once against all remaining cells, reuses that cell-type DEG list in
every signaling direction, and tests every unordered DEG pair within each
receiver before intracellular prior directions are overlaid. Differential
genes are selected using the raw Wilcoxon p-value.

## Workflow

For every ordered sender-to-receiver direction, CellSigN performs:

1. Cell and gene filtering, optional total-count normalization and `log1p`.
2. Selection of up to 5,000 highly variable genes across all cells.
3. One-versus-rest Wilcoxon differential expression over all retained HVGs,
   calculated once for every eligible cell type.
4. Selection of every positive-logFC gene with raw Wilcoxon p-value below
   `--de-threshold`.
   No Top-N gene limit is applied after the HVG step.
5. Construction of every unordered pair among the receiver DEGs.
6. Receiver-cell Pearson correlation and a receiver-specific nested cubic
   model test for every receiver DEG pair. Each edge family is BH corrected.
7. Effect-size calculation and receiver-cell bootstrap stability assessment.
8. Retention of evidence-supported positive-linear, negative-linear, or
   nonlinear edges, followed by assignment of continuous edge costs.
9. Prior-based annotation of receptor and TF roles and known intracellular
   directions. Unannotated association edges are represented in both
   directions during path search; no edge-source class is assigned.
10. Weighted `k` shortest simple directed paths for every eligible R-TF pair.
11. Connection of sender-up ligands to receptors and receiver TFs to
    receiver-up targets. Reversing the ordered pair produces the opposite
    direction independently.

Candidate edges that fail the relation criteria remain in the evidence table
with `Edge_Retained = False`; they do not enter the path graph. Prior knowledge
does not add or remove statistical edges. It supplies node roles and known edge
directions after the undirected receiver association network has been built.

## Statistical definitions

For each eligible cell type, DE is calculated once with a Wilcoxon
rank-sum test comparing that type with all other cells. The full table contains
`P_Value`. A positive-logFC gene is selected when its raw p-value is below
`--de-threshold`.
All qualifying genes are retained; no post-DE Top-N limit is used.

For candidate edge `(u,v)` in receiver `B`, Pearson evidence is

```text
r_B(u,v) = cor(x_u, x_v | cell type = B)
```

and its raw p-values are BH corrected across all receiver DEG pairs.
A linear edge passes when `Q_Corr < alpha` and
`abs(R_Receiver) >= min_abs_r`; its sign determines positive or negative type.

For nonlinear evidence, reduced and full models are compared within receiver
cells in both predictive directions. For `u -> v`, the reduced model is
`v = beta0 + beta1*u`, and the full model adds `u^2` and `u^3`; the reverse
direction is fitted analogously.
The two directional p-values are combined as

```text
p_nonlinear = min(1, 2 * min(p_u_to_v, p_v_to_u))
```

and BH corrected across candidate edges. A nonlinear edge passes when
`Q_Nonlinear < alpha` and `Delta_R2 >= min_delta_r2`. If both linear and
nonlinear criteria pass, the edge is labelled linear.

Relation q-values and effect sizes form a relation score, and bootstrap
recovery frequency forms a stability score. Gene-level DE is used only to
define the receiver gene set and is not included in edge cost. With
user-configurable nonnegative weights:

```text
Data_Cost = 1 - weighted_mean(Relation_Score, Stability)
Final_Cost = Data_Cost + hub_penalty_weight * Hub_Penalty + hop_penalty
Path_Cost = sum(Final_Cost for all intracellular edges in the path)
```

Lower cost indicates stronger combined evidence. `--k-paths` controls how many
lowest-cost simple paths are retained per R-TF pair. There is no maximum path
length. Direct paths between different receptor and TF genes are allowed; a
receptor is excluded as its own zero-length TF endpoint.

Path and receptor permutation p-values use a fixed graph topology and shuffled
edge costs. The pair statistic is the rank-1 (lowest-cost) R-TF path, so all
reported top-k alternatives for that pair share its `Path_P` and `Path_Q`.
The receptor statistic is its best reachable TF. Pair-level `Path_Q` and
receptor-level `Receptor_Q` are reported for diagnosis; paths are not
automatically deleted by those columns.

## Bundled reference data

The intracellular reference is distributed as
`reference_data/human/intracellular_network.txt.gz` because the uncompressed
file exceeds GitHub's 100 MB single-file limit. CellSigN reads the gzip file
directly; manual decompression is unnecessary. See `DATA_PROVENANCE.md` and
`SOURCE_PROVENANCE.md` for sources, limitations, row counts, and checksums.

## Examples and tests

`Example data/` contains a processed 515-cell by 5,000-gene expression matrix
and matching cell-type labels derived from the primary breast cancer dataset
of Chung et al. (2017), available from GEO under accession
[GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688).
The text matrix is an export of the supplied processed h5ad `X` matrix; it is
not the raw GEO matrix. `Example Result/` contains two representative candidate
pathway tables for B-cell-to-malignant and malignant-to-B-cell directions.
These are illustrative outputs rather than experimentally validated pathways.
The input files and their provenance are described in `DATA_PROVENANCE.md`.

`tests/data/` contains a separate, small synthetic fixture used by the
automated end-to-end test. It is not a subset of the Chung dataset.

## Installation

Python 3.9 or newer is required.

```bash
python -m venv .venv
```

Windows CMD:

```cmd
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Linux/macOS:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

## Running both stages

For the default malignant-versus-other bidirectional scope:

```bash
python -m main --stage all \
  --expression /path/data.h5ad \
  --cell-type-column celltype \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --malignant-label Malignant \
  --hvg-top-genes 5000 \
  --de-threshold 0.05 \
  --k-paths 5 \
  --path-permutations 500 \
  --output-dir Results
```

For selected directions, repeat `--pair`:

```bash
python -m main --stage all \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --pair T_cell:Malignant \
  --pair Malignant:T_cell \
  --output-dir Results
```

## Running the two stages separately

The desired `--pair` selection must be supplied in the evidence stage. It
determines the signaling directions. One-vs-rest DE is still calculated once
for every eligible cell type in the dataset. If `--pair` is omitted, the
malignant-other default is saved.

```bash
python -m main --stage evidence \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --pair T_cell:Malignant \
  --output-dir Results

python -m main --stage paths \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --pair T_cell:Malignant \
  --k-paths 5 \
  --path-permutations 500 \
  --output-dir Results
```

The expression filename is appended automatically, for example
`Results_Data_Chung2017_Breast_all`. In a paths-only run, `--expression` only
resolves this output name; the expression matrix is not reread.

## Main parameters

| Parameter | Default | Function |
|---|---:|---|
| `--hvg-top-genes` | 5000 | global HVG count |
| `--de-threshold` | 0.05 | raw Wilcoxon p-value threshold for DE genes |
| `--alpha` | 0.05 | Pearson and nonlinear BH threshold |
| `--min-abs-r` | 0.1 | minimum absolute Pearson effect |
| `--min-delta-r2` | 0.01 | minimum nonlinear incremental R-squared |
| `--stability-resamples` | 20 | receiver-cell bootstrap resamples |
| `--k-paths` | 5 | paths retained per R-TF pair |
| `--path-permutations` | 500 | fixed-topology edge-cost permutations; 0 disables |
| `--relation-weight` | 0.6 | pair relation contribution |
| `--stability-weight` | 0.2 | bootstrap stability contribution |
| `--hub-penalty-weight` | 0.1 | hub penalty multiplier |
| `--hop-penalty` | 0.05 | cost added per intracellular edge |
| `--no-normalize` | off | skip total-count normalization |
| `--no-log1p` | off | skip log transformation |
| `--sample-column` | none | record h5ad sample metadata in the manifest |

## Output structure

```text
Results_<expression>/
|-- <Sender>_to_<Receiver>_pathway.txt
|-- gene_lists/
|   `-- <CellType>_de_genes.txt
|-- evidence/
|   |-- <CellType>_one_vs_rest_de.txt
|   |-- <Receiver>_gene_evidence.txt
|   |-- <Receiver>_edge_evidence.txt
|   |-- <Sender>_to_<Receiver>_rmtf_paths.txt
|   |-- <Sender>_to_<Receiver>_path_edges.txt
|   `-- <Sender>_to_<Receiver>_rmtf_subnetwork_edges.txt
`-- other/
    |-- processed_gene_list.txt
    |-- differential_gene_counts.txt
    |-- test_counts.txt
    `-- run_manifest.json
```

The first five pathway columns remain `Ligand`, `Receptor`, `Mediator`, `TF`,
and `Target`. Commas inside `Ligand` or `Target` denote alternative prior
partners attached to the same intracellular path. Commas inside `Mediator`
denote sequential nodes in path order.

## Reproducibility and interpretation

The manifest records inputs, SHA-256 hashes, parameters, software versions,
cell/sample dimensions, multiple-testing families, and output counts. The
console and `differential_gene_counts.txt` report one DEG count per eligible cell type;
`test_counts.txt` references the same counts for each signaling direction.
Cell types with fewer than five cells are excluded from one-vs-rest tests.

CellSigN outputs candidate association-supported paths. One-vs-rest DE, Pearson
association, nonlinear model comparison, and path permutations do not by
themselves establish causal signaling. Independent biological validation is
required for causal claims.
