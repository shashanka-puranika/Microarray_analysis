# NGS_analysis — RNA-seq pipeline & DE analysis

Overview
--------
This repository contains an end-to-end RNA-seq pipeline (shell) and an edgeR-based differential expression analysis script (R). The pipeline downloads example FASTQ data, performs QC, trimming, alignment (HISAT2), read counting (featureCounts) and generates a counts file. The R script performs normalization, QC plots and differential expression testing using edgeR.

Repository files of interest
---------------------------
- `rnaseq.sh`  
  Shell pipeline that:
  - creates directory layout (rnaseq_analysis/...)
  - downloads example Arabidopsis SRR FASTQs
  - downloads reference genome & annotation (GFF) and converts to GTF (uses `gffread` or downloads pre-converted GTF)
  - runs FastQC / MultiQC, Trimmomatic, HISAT2 index & alignment, samtools sorting/indexing
  - counts reads with `featureCounts` and writes `counts/gene_counts.txt`
  Note: the script creates `rnaseq_analysis` directory and operates inside it.

- `gene_counts_clean.matrix.tsv`  
  A precomputed gene-by-sample counts matrix (rows = genes, columns = sample IDs). Example header: `Geneid  SRR4420293  SRR4420294  SRR4420297  SRR4420298`.

- `design.txt`  
  Simple sample design file used by the R analysis. Must contain (at least) the columns `Sample` and `Group`. Example:
  ```
  Sample  Group
  SRR4420293  control
  SRR4420294  control
  SRR4420297  test
  SRR4420298  test
  ```

- `edgeR_full_analysis_with_normalization.R`  
  edgeR-based DE workflow. Top of the script defines configurable variables:
  - `counts_file` (default: `"gene_counts_clean.matrix.tsv"`)
  - `design_file` (default: `"design.txt"`)
  - `outdir` (plots & HTML index)
  - normalization options (`norm_method`, `use_voom`) and volcano thresholds
  The script:
  - loads counts and design, checks sample name consistency
  - filters low-count genes, normalizes (TMM by default)
  - runs edgeR (glmQLFit / glmQLFTest) with an automatically chosen contrast (works when there are exactly 2 groups)
  - writes `edgeR_results.csv` and produces many QC plots + interactive volcano saved to `outdir`

Requirements
------------
System (command-line tools)
- bash / coreutils
- wget or curl (used in `rnaseq.sh`)
- FastQC
- MultiQC
- Trimmomatic
- HISAT2 (hisat2-build, hisat2)
- samtools
- featureCounts (subread)
- gffread (optional; if not available the script downloads a pre-converted GTF)

R (packages)
The R script will attempt to install missing packages if needed. Required packages include:
- edgeR
- limma
- pheatmap
- ggplot2
- ggrepel
- plotly
- htmlwidgets
- viridis
- dplyr

Quickstart
----------
1. Run the shell RNA-seq pipeline (creates rnaseq_analysis/ and performs alignment & counting)
   - From repository root:
     - Make the script executable and run it:
       ```bash
       chmod +x rnaseq.sh
       ./rnaseq.sh
       ```
     - Or run directly with bash:
       ```bash
       bash rnaseq.sh
       ```
   - Note: the script downloads data and builds HISAT2 indexes; ensure you have enough disk space and required command-line tools installed.

2. Run the R differential expression analysis
   - Confirm `gene_counts_clean.matrix.tsv` and `design.txt` are present in your working directory (or update `counts_file`/`design_file` variables at the top of the R script).
   - From R / RStudio:
     ```r
     # open the script in RStudio and run, or
     source("edgeR_full_analysis_with_normalization.R")
     ```
   - Or run non-interactively:
     ```bash
     Rscript edgeR_full_analysis_with_normalization.R
     ```
   - Results and plots are written to the `outdir` you set in the script (default: `plots_edger`). An `index.html` file is generated to browse plots.

Important notes & troubleshooting
---------------------------------
- Sample name consistency: sample names in `design.txt` must exactly match column names in the counts matrix. The R script checks this and will stop if samples are missing.
- Groups / contrast: the R script automatically creates a contrast when there are exactly two groups. If you have more than two groups, edit the script to set an explicit contrast (see the `makeContrasts()` line).
- GFF -> GTF conversion: `rnaseq.sh` will use `gffread` if available. If not, it downloads a pre-converted GTF from Ensembl plants. Verify that the GTF matches your desired annotation version.
- If you change filenames or paths, update the top variables in `edgeR_full_analysis_with_normalization.R` accordingly.
- The shell pipeline currently creates the analysis folder `rnaseq_analysis` and runs inside it — be mindful of relative paths when running the R script (you can either copy the counts & design into the repo root or point the R script to the files inside `rnaseq_analysis/counts`).

Outputs you can expect
----------------------
- rnaseq pipeline
  - `rnaseq_analysis/` with subfolders: `raw_data`, `qc/`, `trimmed`, `reference`, `alignment`, `counts`, `results`
  - `counts/gene_counts.txt` (featureCounts output)
- edgeR R script
  - `edgeR_results.csv` (DE results)
  - Normalized counts CSVs in `outdir` (e.g., `normalized_CPM.csv`, `voom_E_logCPM.csv` if used)
  - A set of QC and DE plots (PNGs) and an interactive HTML volcano in `outdir`
  - `index.html` in `outdir` for quick browsing of generated plots

Customizing
-----------
- To analyze your own data, replace the FASTQ downloads in `rnaseq.sh` with your own URLs or move your FASTQs to `rnaseq_analysis/raw_data/` and adjust sample names in the script loop.
- To change normalization or other analysis options, edit the variables near the top of `edgeR_full_analysis_with_normalization.R` (e.g., `norm_method`, `use_voom`, volcano thresholds, `top_n_heatmap`).

Author / contact
----------------
Author: Shashanka Puranika K

License
-------
(Choose a license and add here if you wish — currently none specified)

Acknowledgements
----------------
This repository uses many community tools (HISAT2, featureCounts, edgeR, etc.) — please cite the corresponding software when publishing results.
