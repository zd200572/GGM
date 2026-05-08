# GGM

A pipeline for **G**oat **G**ut and **M**icrobiome analysis, comprising three modules: MAG generation, antimicrobial peptide (AMP) identification, and visualization.

## Directory Structure

```
GGM/
├── MAG-Generation/          # Metagenome-assembled genome generation pipeline
├── AMP-Identification/      # Antimicrobial peptide identification pipeline
└── R-script-forFigure/      # R scripts for statistical analysis and visualization
```

---

## MAG-Generation

Pipeline for metagenomic raw read processing, assembly, and binning to recover metagenome-assembled genomes (MAGs).

### Scripts

| File | Description |
|------|-------------|
| `1.qc-assembly.sh` | Read QC, host removal, and assembly pipeline |
| `2.binning.sh` | MAG binning pipeline |
| `remove-host.pl` | Remove host-contaminated reads from paired FASTQ |
| `renamefa.pl` | Rename FASTA headers with prefix and sequential numbering |
| `deal_fa.pl` | Multi-function FASTA utility (20 operations: length/GC stats, rename, filter, reverse complement, N50 calculation, etc.) |

### Workflow

1. **QC & Assembly** (`1.qc-assembly.sh`):
   - Adapter trimming with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   - Host read removal with [BWA](http://bio-bwa.sourceforge.net/) + `remove-host.pl`
   - Assembly with [MEGAHIT](https://github.com/voutcn/megahit)
   - Contig renaming and filtering (>10kb)

2. **Binning** (`2.binning.sh`):
   - Read mapping back to assembly with BWA
   - Depth calculation with [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) JGI tool
   - Binning with MetaBAT2

### Input

- Paired FASTQ files (`1.fastq.gz`, `2.fastq.gz`)
- Adapter file (`adapters.fa`)
- Host genome (`host.fa`)

---

## AMP-Identification

Pipeline for identifying antimicrobial peptides (AMPs) from protein sequences using an ensemble of five prediction tools.

### Scripts

| File | Description |
|------|-------------|
| `AMP-Identify.pl` | Master orchestrator: filters peptides (10-50aa), generates run scripts for all 5 tools |
| `AmPEP-info.pl` | Parse [AmPEP](https://github.com/AbdulhakimDA/AmPEP) prediction output |
| `APIN-info.pl` | Parse [APIN](https://github.com/biomedical-astro/APIN) prediction output |
| `ampir-info.pl` | Parse [ampir](https://github.com/Legana/ampir) prediction output |
| `amplify-info.pl` | Parse [AMPlify](https://github.com/bcgsc/AMPlify) prediction output |
| `ampscan-info.pl` | Parse [AMPscanner-v2](https://github.com/TangYiJing/AMPscanner-v2) prediction output |
| `info-merge.pl` | Merge results from all tools and apply score thresholds |
| `work.sh` | Demo runner script |
| `demo.list` | Example input list (tab-delimited: sample ID + FASTA path) |

### AMP Prediction Tools Used

 AmPEP, APIN, ampir, AMPscanner-v2，AMPlify 

### Input

Tab-delimited list file (see `demo.list`):
```
SampleID    /path/to/protein.faa
```

### Output

Per-sample filtered FASTA, tool-specific prediction results, and merged AMP results (`<sample>.AMP.results`) with columns: ID, Peptide, Method, AMP-score.

---

## R-script-forFigure

R scripts for statistical analysis and publication-quality visualization of diversity, abundance, and correlation data.

### Scripts

| File | Input | Output | Description |
|------|-------|--------|-------------|
| `alpha.boxplot.r` | `diversity.xls`, `group.xls` | `Diff_alpha.pdf` | Alpha diversity boxplot (Richness/Shannon/Simpson) with statistical annotations |
| `diff_alpha.r` | — | — | `diff.alpha()` function: t-test (2 groups) or ANOVA+Tukey HSD (>2 groups) for alpha diversity |
| `barplot.top10.r` | `phylum.abundance.xls`, `group.xls` | `phylum_top10.pdf` | Stacked alluvial barplot of top-10 phyla |
| `heatmap.50.r` | `rpkm.ampID.xls`, `group.xls` | `heatmap.pdf`, `heatmap.xls` | Z-score heatmap of top-50 AMPs with significance letters |
| `pca.r` | `rpkm.xls`, `pca.group.xls` | `pca.pdf` | PCA biplot (PC1 vs PC2) with 95% confidence ellipses |
| `correlation.r` | `correlation.ampID.xls`, `correlation.meta.xls` | `top50_correlation.pdf`, correlation tables | Spearman correlation dotplot between top-50 AMPs and metabolites |
| `correlation.dotplot.r` | — | — | `correlation.dotplot()` function for generating correlation dotplots |

### Required R Packages

`ggplot2`, `vegan`, `reshape2`, `tidyverse`, `dplyr`, `multcomp`, `ggalluvial`, `ggrepel`, `pheatmap`, `patchwork`, `linkET`

---

## Citation

If you use this pipeline, please cite the relevant tools listed above.
