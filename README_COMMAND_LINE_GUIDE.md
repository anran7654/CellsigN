# CellSigN Command-Line Guide

See [README.md](README.md) for the statistical workflow and output definitions.

## One-stage run

```bash
python -m main --stage all \
  --expression /path/data.h5ad \
  --cell-type-column celltype \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --malignant-label Malignant \
  --hvg-top-genes 5000 \
  --de-threshold 0.05 \
  --alpha 0.05 \
  --min-abs-r 0.1 \
  --min-delta-r2 0.01 \
  --stability-resamples 20 \
  --k-paths 5 \
  --path-permutations 500 \
  --output-dir Results \
  --seed 42
```

`--de-threshold` is the raw Wilcoxon p-value cutoff for positive-logFC genes.

## Selected ordered pairs

```bash
python -m main --stage all \
  --expression /path/data.h5ad \
  --intracellular-prior reference_data/human/intracellular_network.txt.gz \
  --ligand-receptor-prior reference_data/human/ligand_receptor.txt \
  --pair T_cell:Malignant \
  --pair Malignant:T_cell \
  --output-dir Results
```

## Separate evidence and path stages

The ordered pair selection belongs to the evidence definition. Use
the same pair in both commands, or omit `--pair` in the paths command to run
all directions stored by the evidence stage.

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

## Windows CMD single-line example

```cmd
".venv\Scripts\python.exe" -m main --stage all --expression "E:\data\data.h5ad" --cell-type-column celltype --intracellular-prior "reference_data\human\intracellular_network.txt.gz" --ligand-receptor-prior "reference_data\human\ligand_receptor.txt" --malignant-label Malignant --hvg-top-genes 5000 --de-threshold 0.05 --alpha 0.05 --k-paths 5 --path-permutations 500 --output-dir Results --seed 42
```

Use `python -m main --help` for the complete current parameter list.
