# 🧬 transcriptomics-rnaseq-pipeline

> Reproducible RNA-seq pipeline from raw FASTQ to pathway enrichment and network analysis — designed to run on Google Colab.

![Platform](https://img.shields.io/badge/platform-Google%20Colab-F9AB00?logo=googlecolab)
![Conda](https://img.shields.io/badge/conda-mamba-green)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)
![Status](https://img.shields.io/badge/status-in%20development-yellow)

---

## Overview

This pipeline processes Illumina paired-end RNA-seq data from SRA through quality control, trimming, alignment, quantification, differential expression analysis, enrichment, and network analysis.

**Dataset:** GSE183117 — Mouse retina P10, rd1 mutation vs wild-type  
**Organism:** *Mus musculus* (GRCm39 / Ensembl 110)  
**Contrast:** rd1 vs WT

```
SRA FASTQ
    │
    ├── [Step 0]  Setup — condacolab + mamba + R packages
    ├── [Step 1]  Download — sra-tools (prefetch + fasterq-dump)
    ├── [Step 2]  QC — FastQC + MultiQC
    ├── [Step 3]  Trimming — fastp (Q≥20, auto adapters)
    ├── [Step 4]  Alignment — HISAT2 (splice-aware, ~8 GB RAM)
    ├── [Step 5]  Quantification — featureCounts (raw integer counts)
    ├── [Step 6]  DESeq2 QC — PCA + Pearson correlation heatmap (all samples)
    ├── [Step 7]  Outlier removal + DESeq2 re-fit
    ├── [Step 8]  Differential expression — DESeq2 (BH/FDR correction)
    ├── [Step 9]  ID conversion — ENSEMBL → symbol + Entrez (org.Mm.eg.db)
    ├── [Step 10] Visualization — plotMA, EnhancedVolcano, DEG heatmap
    ├── [Step 11] Enrichment — enrichGO + gost + GSEA (dotplot + table each)
    └── [Step 12] Network — STRINGdb + igraph (hub genes + PPI network)
```

---

## Repository Structure

```
transcriptomics-rnaseq-pipeline/
├── transcriptomics_pipeline_v2.ipynb   # Main notebook (Google Colab)
├── README.md
├── .gitignore
└── LICENSE
```

---

## Tool choices and rationale

| Step | Tool chosen | Alternative | Reason |
|------|------------|-------------|--------|
| Trimming | **fastp** | Trimmomatic | 2–4× faster, auto adapter detection, built-in QC report |
| Alignment | **HISAT2** | STAR | ~4× less memory (8 GB vs 30 GB), equally accurate for known splicing |
| Quantification | **featureCounts** | HTSeq | 10–20× faster, multi-threaded, single command for all samples |
| Normalisation | **DESeq2 median-of-ratios** | manual RPM/FPKM | Correct for count data; never pre-normalise before DESeq2 |
| p-value correction | **BH (FDR)** | Bonferroni | Controls false discovery rate; standard for RNA-seq |
| Enrichment A | **enrichGO** (up/down) | — | ORA with controlled background; reproducible via Bioconductor |
| Enrichment B | **gost** (up/down) | — | Multi-database in one call; auto-updated via g:Profiler API |
| Enrichment C | **GSEA** (all genes ranked) | — | No cutoff needed; captures subtle coordinated signals |
| Network | **STRINGdb + igraph** | — | Identifies hub genes for experimental validation |

---

## Requirements

- **Google Colab** (recommended — free tier sufficient, ~12 GB RAM)
- Or: Linux / WSL2 with Mamba installed

### Tools installed automatically

| Tool | Conda environment |
|------|------------------|
| sra-tools | `sra` |
| FastQC, MultiQC, fastp | `qc` |
| HISAT2, Samtools, featureCounts (subread) | `hisat2` |

### R packages installed automatically

`DESeq2`, `EnhancedVolcano`, `clusterProfiler`, `enrichplot`, `org.Mm.eg.db`, `AnnotationDbi`, `ReactomePA`, `STRINGdb`, `ggplot2`, `pheatmap`, `gplots`, `ggrepel`, `gprofiler2`, `igraph`, `ggraph`, `writexl`, `dplyr`, `RColorBrewer`

---

## Input samples

| SRR | Condition | Replicate |
|-----|-----------|-----------|
| SRR15679348 | WT | 1 |
| SRR15679349 | WT | 2 |
| SRR15679350 | WT | 3 |
| SRR15679351 | WT | 4 |
| SRR15679346 | rd1 | 1 |
| SRR15679347 | rd1 | 2 |

All files are downloaded automatically inside the notebook from NCBI SRA.

---

## Usage

### Google Colab (recommended)

1. Open [Google Colab](https://colab.research.google.com/)
2. Upload `transcriptomics_pipeline_v2.ipynb`
3. Run **Step 0.1** and wait for the kernel to restart
4. Continue running cells sequentially

### Local (Linux / WSL2)

Replace `%%bash` cells with terminal commands and `mamba run -n <env>` with `conda activate <env>`.

---

## Output structure

```
figures/
├── qc/           # FastQC, MultiQC (pre/post-trim), fastp HTML, HISAT2 summaries
├── deseq2/       # 01_pca_all · 02_heatmap_all · 03_pca_filt · 04_heatmap_filt
│                 # 05_plotMA · 06_EnhancedVolcano · 07_heatmap_DEGs
├── enrichment/   # 11A_enrichGO · 11B_gost · 11C_GSEA (dotplots)
└── network/      # 12_PPI_network

tables/
├── DESeq2_all_results_rd1.csv        # All tested genes (Mutation = RD1)
├── 11A_enrichGO_results.xlsx         # enrichGO up + down (Mutation = RD1)
├── 11B_gost_results.xlsx             # gost up + down (Mutation = RD1)
├── 11C_GSEA_GO_results.xlsx          # GSEA GO results (Mutation = RD1)
└── 12_hub_genes.xlsx                 # Top 15 hub genes (Mutation = RD1)
```

---

## Key analysis decisions

### Why raw counts in DESeq2?
DESeq2 requires **raw integer counts** as input and performs normalisation internally via the median-of-ratios method. Pre-normalised values (RPKM, FPKM, TPM) or rounded decimals produce incorrect p-values and fold changes.

### Why re-run DESeq2 after outlier removal?
Size factors and dispersion estimates depend on all samples. Removing an outlier changes both — the model must be re-fitted on the clean dataset.

### Why use `padj` (not `pvalue`) in the volcano plot?
The raw `pvalue` is not corrected for multiple testing. With ~20,000 genes tested, hundreds will appear significant by chance. `padj` (BH/FDR) is the correct metric for significance calls and all downstream analyses.

### Why separate up/down in ORA (enrichGO, gost)?
ORA tests whether a gene set is over-represented in a list. Mixing up- and down-regulated genes in the same list can cancel out real signals and produces biologically ambiguous results. Separating them allows identifying which pathways are **activated** vs **repressed**.

### GSEA ranking metric
`sign(log2FC) × −log10(padj + ε)` captures both the direction of regulation and statistical confidence. This is preferred over log2FC alone for GSEA.

---

## Tools and References

| Tool | Reference |
|------|-----------|
| fastp | [Chen et al., 2018](https://doi.org/10.1093/bioinformatics/bty560) |
| HISAT2 | [Kim et al., 2019](https://doi.org/10.1038/s41587-019-0201-4) |
| featureCounts | [Liao et al., 2014](https://doi.org/10.1093/bioinformatics/btt656) |
| DESeq2 | [Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8) |
| EnhancedVolcano | [Blighe et al.](https://github.com/kevinblighe/EnhancedVolcano) |
| clusterProfiler | [Wu et al., 2021](https://doi.org/10.1016/j.xinn.2021.100141) |
| gprofiler2 | [Kolberg et al., 2021](https://doi.org/10.1093/nar/gkab301) |
| STRINGdb | [Szklarczyk et al., 2021](https://doi.org/10.1093/nar/gkaa1074) |
| org.Mm.eg.db | [Bioconductor](https://bioconductor.org/packages/org.Mm.eg.db) |

---

## Notes

- Run **Step 0.1** first and wait for the kernel to restart before continuing.
- **Strandedness** (`-s` in featureCounts, `--rna-strandness` in HISAT2): use `-s 2` / `RF` for TruSeq Stranded kits. If unsure, run RSeQC `infer_experiment.py` on one BAM file first.
- The HISAT2 index download (~8 GB) is the longest step (~15 min on Colab).
- All result tables include a `Mutation = "RD1"` column for downstream merging.
- This pipeline was developed and tested on Google Colab free tier.

---

## Author

**Gabriel Ibovi**  
Bioinformatics | Transcriptomics | Multi-omics Data Analysis  
🔗 [github.com/gabrielibovi-bioinfo](https://github.com/gabrielibovi-bioinfo)

---

## License

MIT License — see the [LICENSE](LICENSE) file for details.
