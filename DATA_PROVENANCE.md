# Reference data provenance

The two tab-delimited reference tables in `reference_data/human` were migrated
without row-level changes from the supplied `CellSign_0413.zip`. The intracellular
table is stored as `.txt.gz` to satisfy GitHub's file-size limit; CellSigN reads
it directly and decompression reproduces the original byte stream. The
ligand-receptor table remains an uncompressed `.txt` file.

The intracellular aggregate contains the columns `Gene1`, `Gene2`, and `Type`. The ligand-receptor aggregate contains `Ligand`, `Receptor`, and `Score`, but this release treats the supplied file as already curated: it reads every nonempty ligand-receptor record from the first two columns and ignores `Score`. Neither legacy file records a source database for every row. Consequently, the repository cannot independently assign an edge-level meaning to the score or verify a claim that every interaction is supported by a particular number of databases.

According to the manuscript supplied with the code, the TF-target component was assembled from TRRUST, TFLink, hTFtarget, HTRIdb, and RegNetwork; other intracellular interactions were obtained through Pathway Commons. The ligand-receptor aggregate was assembled from CellPhoneDB, CellTalkDB, OmniPath, and the Lewis Lab collection. Users should consult and comply with each upstream resource's current terms before redistributing derived data.

The GPL-3.0 license in this repository applies to the source code. It does not replace the terms of the upstream databases.

## Real-data example

The primary breast cancer single-cell dataset of Chung et al. (2017) is
publicly listed by GEO as accession
[GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688).
GEO describes 549 single-cell profiles before quality filtering and 515
retained single-cell profiles. The repository's `Example data/` files were
exported from the processed `Data_Chung2017_Breast_all.h5ad` used for the
CellSigN analysis, retaining its 515 cells, 5,000 selected genes, cell IDs,
cell-type labels, and `X` expression values. They are not the original GEO
raw/processed-download files. The original h5ad is not included in the
repository.

`Example Result/` contains representative B-cell-to-malignant and
malignant-to-B-cell pathway tables from an existing CellSigN run on that
processed h5ad. The example outputs were selected from a run with raw
one-versus-rest Wilcoxon p-value threshold 0.05, `k_paths=1`, 500 path
permutations, seed 42, and normalization/log1p disabled. The delimited
expression export preserves the h5ad values as round-trippable text; ordinary
CSV parsers may introduce floating-point differences at approximately
machine precision. These tables are candidate associations, not experimental
validation.

The small synthetic matrices under `tests/data/` are test fixtures and have
no biological provenance claim.
