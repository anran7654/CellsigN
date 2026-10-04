# CellSigN 0.3.5

CellSigN infers candidate bidirectional `Ligand-Receptor-Mediator-TF-Target`
paths from single-cell RNA-seq data. Version 0.3.5 tests each eligible
cell type once against all remaining cells and reuses that cell-type gene list
in every signaling direction. Differential genes are selected using the raw
Wilcoxon p-value. When fewer than the requested number pass, CellSigN can
supplement the downstream list with the lowest-p-value non-DE genes.

## Workflow

For every ordered sender-to-receiver direction, CellSigN performs:

1. Cell and gene filtering, optional total-count normalization and `log1p`.
2. Selection of up to 5,000 highly variable genes across all cells.
3. One-versus-rest Wilcoxon differential expression over all retained HVGs,
   calculated once for every eligible cell type.
4. Selection of every positive-logFC gene with raw Wilcoxon p-value below
   `--de-threshold`. If fewer than `--de-gene-target` genes pass, genes with
   p-values above the threshold are added in ascending p-value order until the
   target is reached. A larger genuine DEG set is never truncated.
5. Construction of every unordered pair among the receiver selected genes.
6. Receiver-cell Pearson correlation and a receiver-specific nested cubic
   model test for every receiver selected-gene pair. Each edge family is BH
   corrected.
7. Effect-size calculation and receiver-cell bootstrap stability assessment.
8. Retention of evidence-supported positive-linear, negative-linear, or
   nonlinear edges, followed by assignment of continuous edge costs.
9. Intersection of evidence-supported edges with the intracellular signaling
   reference. Only intersected edges enter path search and follow reference
   directions; receptor and TF roles also come from the prior.
10. Weighted `k` shortest simple directed paths for every eligible R-TF pair.
11. Connection of sender-up ligands to receptors and receiver TFs to
    receiver-up targets. Reversing the ordered pair produces the opposite
    direction independently.

Candidate edges that fail the relation criteria remain in the evidence table
with `Edge_Retained = False`. Evidence-supported edges absent from the
intracellular reference also remain in that table. Neither group enters the
path graph, which requires both statistical support and a signaling-prior
annotation. The prior supplies node roles and directions for those edges.

## Statistical definitions

For each eligible cell type, DE is calculated once with a Wilcoxon
rank-sum test comparing that type with all other cells. The full table contains
`P_Value`. A positive-logFC gene is selected when its raw p-value is below
`--de-threshold`. All qualifying genes are retained. If this DEG set contains
fewer than `--de-gene-target` genes (default 500), CellSigN supplements it with
genes having `P_Value > --de-threshold`, ranked by ascending raw p-value, until
the target is reached. The target is a supplementation ceiling rather than an
upper bound: for example, 614 DE genes remain 614, whereas 25 DE genes are
supplemented to 500 when enough eligible genes exist. If the entire tested gene
universe contains at most the target number, all tested genes are retained.
Set `--de-gene-target 0` to disable supplementation.

For candidate edge `(u,v)` in receiver `B`, Pearson evidence is

```text
r_B(u,v) = cor(x_u, x_v | cell type = B)
```

and its raw p-values are BH corrected across all receiver selected-gene pairs.
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
python -m cellsign --stage all \
  --expression /path/data.h5ad \
  --cell-type-column celltype \
  --intracellular-prior reference_data/human/intracellular_network.txt \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --malignant-label Malignant \
  --hvg-top-genes 5000 \
  --de-threshold 0.05 \
  --de-gene-target 500 \
  --k-paths 5 \
  --path-permutations 500 \
  --output-dir Results
```

For selected directions, repeat `--pair`:

```bash
python -m cellsign --stage all \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt \
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
python -m cellsign --stage evidence \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --pair T_cell:Malignant \
  --output-dir Results

python -m cellsign --stage paths \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt \
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
| `--de-gene-target` | 500 | supplement smaller DEG sets to this downstream target; 0 disables |
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
