# Tissue-Specific Cancer Risk: Chromosomal Deletion Analysis

Analysis code accompanying the study of **tissue-specific patterns of chromosomal arm-level copy number deletions** (13q and 17q) in *BRCA1/BRCA2* germline carriers, integrating bulk tumour cohorts (TCGA, ICGC) and single-cell whole-genome sequencing (scWGS) across benign, pre-malignant, and malignant breast and ovarian tissues.

---

## Table of Contents

- [Background](#background)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Workflow Overview](#workflow-overview)
- [Scripts Reference](#scripts-reference)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [License](#license)

---

## Background

Germline *BRCA1* and *BRCA2* mutations predispose carriers to breast and ovarian cancer, in part through somatic loss of the wild-type allele on chromosomes 17q and 13q respectively. This project characterises:

1. **Cohort-level prevalence** — the frequency of 13q/17q arm-level deletions across cancer types using TCGA and ICGC bulk-tumour data, after excluding whole-genome duplication (WGD) events.
2. **Single-cell clonal evolution** — haplotype-specific copy number profiling in individual *BRCA1/BRCA2* carrier samples spanning normal breast epithelium → atypical lobular hyperplasia (ALH) → lobular/ductal carcinoma in situ (LCIS/DCIS) → invasive carcinoma.
3. **Tissue specificity** — whether deletion enrichment is cancer-type-specific using permutation-based statistical testing.

---

## Repository Structure

```text
TissueSpecificCancerRisk/
├── src/                        # All analysis scripts (R and Python)
│   ├── 0.Process_ICGC_toNoWGD.R
│   ├── 0.Process_TCGA_toNoWGD.R
│   ├── 1.chr_arm_del_permutation.py
│   ├── 2.Prevalence_plot_desc.R
│   ├── aneuploidy_analysis.R
│   ├── clustering.R
│   ├── heatmap_plot.R
│   ├── make_heatmap.R
│   ├── util.R
│   ├── col_palettes.R
│   └── ...                     # see Scripts Reference below
├── jobs/                       # HPC job submission scripts (LSF/bsub)
├── README.md
└── .gitignore
```

> Raw data, large intermediate files, results, and figure outputs are not tracked in this repository (see [Data Availability](#data-availability)).

---

## Dependencies

### R (>= 4.2)

| Package | Purpose |
| --- | --- |
| `tidyverse` | Data wrangling |
| `data.table` / `vroom` | Fast reading of large CNV files |
| `ComplexHeatmap` | Copy-number heatmaps with tree sidebars |
| `ggplot2` / `ggnewscale` | Visualization |
| `ape` / `ggtree` | Phylogenetic tree I/O and plotting |
| `lme4` / `sandwich` / `lmtest` | Mixed-effects models, robust inference |
| `dbscan` | HDBSCAN clustering |
| `umap` | UMAP dimensionality reduction |
| `matrixStats` | Efficient matrix operations |
| `here` | Portable project-relative paths |
| `openxlsx` / `readxl` | Excel file I/O |

Install all R dependencies:

```r
pkgs <- c(
  "tidyverse", "data.table", "vroom",
  "ggplot2", "ggnewscale", "ape",
  "lme4", "sandwich", "lmtest",
  "dbscan", "umap", "matrixStats",
  "here", "openxlsx", "readxl"
)
install.packages(pkgs)
BiocManager::install(c("ComplexHeatmap", "ggtree"))
```

### Python (>= 3.9)

```text
pandas
numpy
statsmodels
scipy
```

```bash
pip install pandas numpy statsmodels scipy
```

---

## Workflow Overview

```text
Raw TCGA / ICGC CNV data
        │
        ▼
0.Process_TCGA_toNoWGD.R          # Filter to primary tumours, remove WGD
0.Process_ICGC_toNoWGD.R
        │
        ▼
1.chr_arm_del_permutation.py      # 10,000 permutations × cancer type × chr arm
                                  # → enrichment p-values + BH-FDR correction
        │
        ▼
2.Prevalence_plot_desc.R          # Prevalence bar plots and summary figures

Single-cell scWGS CNV data (per sample)
        │
        ▼
cnv_processor.R / read_cnv_data.R # Wide-format conversion, per-patient filtering
        │
        ├── clustering.R           # UMAP + HDBSCAN → clone assignments
        ├── deletion-threshold.Rmd # Threshold calibration for deletion calls
        ├── brca-del-17q13q.Rmd   # Per-sample 13q/17q deletion analysis
        │
        ▼
make_heatmap.R / heatmap_plot.R   # Publication heatmaps with phylogenetic trees
        │
        ▼
clonality-analysis.Rmd            # Integrated clonal structure analysis
aneuploidy_analysis.R             # Aneuploidy statistics + mixed-effects models
```

---

## Scripts Reference

### Data Preprocessing

| Script | Description |
| --- | --- |
| `0.Process_TCGA_toNoWGD.R` | Filter TCGA CNV data to primary tumours without WGD; retain 13q/17q deletion calls |
| `0.Process_ICGC_toNoWGD.R` | Same preprocessing for ICGC dataset |
| `cnv_processor.R` | Convert long-format CNV bins to wide sample × bin matrix; supports `--patient` / `--sample_id` filters |
| `read_cnv_data.R` | Utility reader mirroring `cnv_processor.R` for in-memory use |

### Statistical Analysis

| Script | Description |
| --- | --- |
| `1.chr_arm_del_permutation.py` | Permutation test (n = 10,000) for chromosome arm deletion enrichment per cancer type; outputs p-values with BH-FDR correction |
| `aneuploidy_analysis.R` | Mixed-effects models (lme4) for aneuploidy burden; robust standard errors via sandwich/lmtest |

### Single-Cell / Clonal Analysis

| Script | Description |
| --- | --- |
| `clustering.R` | UMAP + HDBSCAN clustering pipeline for copy number profiles (supports haplotype-specific CN) |
| `cluster_A_B.R` | Allele-specific (A/B) copy number clustering |
| `clustering_cells.R` | Cell-level clustering utilities |
| `cntob.R` | Copy number to binary deletion/gain conversion |
| `deletion-threshold.Rmd` | Notebook for calibrating deletion-call thresholds |
| `brca-del-17q13q.Rmd` | Per-sample analysis of 13q/17q deletions in *BRCA1/2* carriers across tissue states |
| `luminal.Rmd` | Luminal breast epithelial cell-specific analysis |
| `clonality-analysis.Rmd` | Integrated clonality analysis across tissue type, chemotherapy history, and tumour stage |
| `data_preparation.Rmd` | Data preparation for SITKA phylogenetic tree inference |

### Visualization

| Script | Description |
| --- | --- |
| `make_heatmap.R` | Driver for publication-quality CNV heatmaps with phylogenetic tree sidebars |
| `heatmap_plot.R` | Core heatmap-building functions (ComplexHeatmap-based) |
| `vizualize_tree.R` | Phylogenetic tree visualization utilities |
| `plot_tree_deletions.R` | Overlay deletion calls onto inferred phylogenetic trees |
| `col_palettes.R` | Colour palettes for CN states (0-11+), HSCN, and deletion overlays |
| `global_aes_out.R` | Global ggplot2 theme and aesthetic defaults |
| `2.Prevalence_plot_desc.R` | Descriptive prevalence plots across cancer types |

### Utilities

| Script | Description |
| --- | --- |
| `util.R` | Core utilities: `createCNmatrix()`, `Mode()`, genome-binning helpers |
| `print.methods.R` | Custom S3 print methods |

---

## Data Availability

Large raw and intermediate data files are not included in this repository. Required data can be obtained from:

- **TCGA** copy number segment data: [GDC Data Portal](https://portal.gdc.cancer.gov/)
- **ICGC** somatic CNV data: [ICGC Data Portal](https://dcc.icgc.org/)
- **Single-cell WGS data**: accession numbers and processed CNV matrices will be deposited in a public repository upon publication (e.g., EGA / GEO)

Once downloaded, place raw files in `data/` following the path conventions used in the preprocessing scripts.

---

## Citation

> *Manuscript in preparation.* Citation will be added upon publication.

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

For questions about the analysis code, please open a [GitHub Issue](../../issues).
