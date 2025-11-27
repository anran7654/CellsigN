# CellsigN
CellSigN, a computational framework for reconstructing complete “ligand–receptor–mediator–TF–target” signaling pathways. CellSigN integrates ligand–receptor interactions, transcriptional regulation, and six types of intracellular signal transduction.
![git push](workflow.png)

As shown in figure, CellSigN begins by constructing cell type–specific gene correlation networks from scRNA-seq data and corresponding cell annotations. It then applies a contrastive-cutting strategy to integrate known TF–target (TF–TG) and receptor–TF interactions, curated from multiple databases, with these correlation networks to generate intracellular signaling subnetworks. By incorporating five major categories of molecular interactions—physical association, complex formation, catalytic regulation, state modulation, and phosphorylation—CellSigN identifies the molecular roles within each subnetwork and classifies them as receptors, intermediate mediators, or TF. The shortest receptor–TF paths are subsequently determined using a breadth-first search (BFS) algorithm, enabling the connection of ligands in sender cells to TF in receiver cells through receptor–mediator–TF cascades. Known ligand–receptor pairs are then mapped onto the cell type–specific networks. Integrating both intracellular and intercellular components, CellSigN ultimately reconstructs complete cross-cell signaling pathways that encompass the ligand–receptor–mediator–TF–target continuum.
### Preparation

CellSigN is a Python script tool. Python environment needs:

- Python 3.12 or above  
- PyTorch 1.4.0 or above  
- numpy 1.26.4  
- pandas 2.2.2  
- scanpy 1.11.0

# Input
The input type of each single-cell sequencing data expression matrix is `​​.csv` or `.txt`, where rows represent cell names and columns represent features (genes).It is important that cellsign also inputs the Meta data file in order to obtain the cell type of the cell in the analysis. The example expression matrix and meta data is available for download under `Example Data/`.
- 1. When running `run.py` from the command line, provide the single-cell sequencing expression matrix path via `--file_path` and the metadata path via `--cell_path`.

- Excerpt from the example matrix：

| Gene      | BC01_02      | BC01_03      | BC01_04      | BC01_05      | BC01_06      |
|-----------|--------------|--------------|--------------|--------------|--------------|
| TSPAN6    | 0            | 0            | 0            | 0            | 0            |
| DPM1      | 3.946252186  | 4.798721645  | 2.07616304   | 5.202344561  | 3.50040145   |
| SCYL3     | 1.454120482  | 0.592937371  | 4.535610877  | 2.006141851  | 4.10844772   |
| C1orf112  | 0            | 2.667806009  | 0.292808481  | 0.230937983  | 0            |
| FGR       | 0            | 0            | 0            | 0            | 0            |
| CFH       | 0            | 0            | 0            | 0            | 0.285478504  |
| FUCA2     | 2.364818474  | 3.619753408  | 0            | 0            | 2.243913826  |


- Meta data：

| cell_name | cell_type |
|-----------|-----------|
| BC01_02   | Malignant |
| BC01_03   | Malignant |
| BC01_04   | Malignant |
| BC01_05   | Malignant |
| BC01_34   | Malignant |
| BC01_50   | Stromal   |
| BC01_53   | Malignant |

- 2.Modify the complete intracellular signaling library file path (`pathway_file`) and the curated ligand–receptor library（`ligand_file`） file path in `run.py`.The librarys are available for download under `Interaction Database/`.

## Parameter Adjustment

### Required parameters

The following arguments must be provided when running `run.py`:

- `--file_path` (str)  
  Path to the gene expression matrix file (e.g. `.txt`).

- `--cell_path` (str)  
  Path to the cell-type annotation file. (e.g. `.txt` or `.csv`).

- `--pathway_file` (str)  
  Path to the pathway interaction file (e.g. `Pathway_Commons.txt`).

- `--ligand_file` (str)  
  Path to the ligand–receptor interaction file (e.g. `ligand-receptor.txt`).

- `--output_dir` (str)  
  Directory for saving all output files. It will be created automatically if it does not exist.

---

### Optional parameters
Detailed parameter descriptions and default values ​​are as follows:
- `--min_cell` (float, default: `0.01`)  
  parameter for gene filtering

- `--min_gene` (float, default: `0.01`)  
  parameter for cell filtering

- `--normalize` (bool, default: `True`)  
  normalize cells

- `--log_trans` (bool, default: `True`)  
  logarithm expression

- `--cell_top_gene` (int, default: `500`)  
  take the top X genes expressed in each cell (default: 1000)

- `--alpha` (float, default: `0.01`)  
  Significant level. Larger alpha leads to more edges (0–1 range)

- `--max_length` (int, default: `10`)  
  Maximum length of the receptor–TF pathway

- `--min_cells` (float, default: `3`)  
  parameter for Filter out cell types larger than a certain number of cells

# Optional usage
- `gene_list =`（This usage is enabled by default）
  You can specify the output path to output the preprocessed gene list

# Contact
For any questions please contact Miss. Yu Zhao (Email: zhaoyu@ncpsb.org.cn).
