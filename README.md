# Cell-Free Multi-Omics Analysis in Gastrointestinal Cancer

<p align="center">
  <img src="docs/images/upstream/CRC_STAD_NC_MultiOmics_Graphical_Abstract.svg" alt="CRC-STAD-NC Multi-Omics Workflow" width="100%">
</p>


A production-grade Nextflow DSL2 pipeline reproducing the computational analyses from [Tao et al. 2023 (*Cell Reports Medicine* 4:101281)](https://doi.org/10.1016/j.xcrm.2023.101281) across four public cell-free omics layers in colorectal cancer (CRC), gastric cancer (STAD), and healthy donors.

---
## Research Team

**Horizon Science Communication LLC — Bioinformatics Team**

- **Hesham Abdullah, PhD**
- **Noura Khaled**

---

## Project Background

Gastrointestinal (GI) cancers, including colorectal cancer (CRC) and stomach adenocarcinoma (STAD), represent major challenges in cancer detection and molecular characterization. Identifying reliable molecular features in blood samples may support the development of less invasive approaches for cancer research and biomarker discovery.

Cell-free DNA (cfDNA) and cell-free RNA (cfRNA) circulating in blood can carry molecular information associated with underlying biological and disease processes. However, individual molecular layers capture different aspects of cancer biology and may provide complementary information.

This project investigates a cell-free multi-omics approach by analyzing publicly available sequencing datasets across colorectal cancer, stomach adenocarcinoma, and healthy donor groups. The study incorporates whole-genome sequencing of cfDNA (cfWGS), methylated cfDNA profiling (cfMeDIP-seq), total cfRNA sequencing, and small cfRNA sequencing for miRNA profiling.

By characterizing these molecular layers, the project aims to investigate cancer-associated genomic, epigenetic, and transcriptomic patterns and explore their potential biological relevance.

---

## Project Objectives

The main objective of this project is to characterize cell-free molecular profiles across colorectal cancer, stomach adenocarcinoma, and healthy donor groups using a multi-omics bioinformatics workflow.

The specific objectives are to:

- Characterize cfDNA genomic alterations, with a focus on copy-number variation (CNV) features.
- Investigate complementary cfDNA fragmentomic features, including end-motif patterns and fragment-size profiles.
- Profile DNA methylation patterns using cfMeDIP-seq data.
- Quantify gene expression from total cfRNA-seq data and investigate relevant RNA features.
- Characterize small cfRNA-derived miRNA profiles and explore their biological relevance.
- Integrate selected genomic, methylation, and gene-expression features to investigate cross-omics relationships and biologically relevant pathways.
- Explore molecular signatures that may support future biomarker discovery in gastrointestinal cancers.
  
---

## Why Multi-Omics Integration?

Cancer involves interconnected genomic, epigenetic, and transcriptomic alterations. Analyzing these molecular layers together may provide a more comprehensive view of cancer biology than examining each layer independently.

This project integrates selected features from three core layers:
- **cfWGS:** Copy-number variations (CNVs).
- **cfMeDIP-seq:** DNA methylation patterns.
- **Total cfRNA-seq:** Gene-expression profiles.

Combining these complementary features may reveal cross-omics relationships and biologically relevant cancer-associated patterns. Fragmentomic and miRNA features are investigated separately to provide additional molecular insights.

---

## Study Questions

- Can liquid biopsy-derived molecular features support the early detection of CRC, particularly at stages I–II, compared with healthy donors?
- Which molecular features distinguish CRC and STAD from healthy donors, and which are cancer-type-specific?
- How are genomic, epigenetic, and transcriptomic alterations related across the three core omics layers?
- Which biological pathways and molecular signatures emerge from multi-omics integration?
- What additional disease-associated information can be captured through cfDNA fragmentomic and small cfRNA-derived miRNA profiles?
  
---

## Study Design

| Parameter | Details |
|---|---|
| **Participants** | 161 individuals: 44 CRC · 36 STAD · 81 healthy donors |
| **Datasets** | 360 published (352 after QC) |
| **GEO Series** | [GSE186607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186607) |
| **Disease focus** | Colorectal cancer (CRC) and gastric cancer (STAD) |

---

## Omics Layers

| Layer | Description | Datasets | BioProject |
|---|---|---|---|
| **cfWGS** | Whole-genome sequencing of cfDNA — copy-number and fragmentomic features | 86 | [PRJNA774626](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA774626) |
| **cfMeDIP-seq** | Antibody-enriched methylated cfDNA — promoter and CpG-island methylation | 98 | [PRJNA774625](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA774625) |
| **Total cfRNA-seq** | Long and short cfRNA — gene expression, splicing, chimeric RNA | 127 | [PRJNA774661](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA774661) |
| **Small cfRNA-seq** | Small RNA / miRNA profiling (QIAseq, UMI-based) | 49 | [PRJNA774739](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA774739) |

---

## Pipeline Architecture

The workflow is implemented in **Nextflow DSL2** and runs on **SLURM HPC**. All four omics branches are executed in parallel within a single `nextflow run` invocation.

```
main.nf
 ├── cfWGS branch          fasterq-dump → FastQC → Trim Galore → Bowtie2 → Samtools → GATK MarkDuplicates → Picard insert size
 ├── cfMeDIP-seq branch    fasterq-dump → FastQC → fastp → Bowtie2 → Samtools → GATK MarkDuplicates → 300-bp bin counts
 ├── Total cfRNA branch    fasterq-dump → FastQC → Cutadapt → STAR (GRCh38) → featureCounts · ERCC (Bowtie2) · rRNA (Bowtie2) → spikein QC
 └── Small cfRNA branch    fasterq-dump → FastQC → Trim Galore → Bowtie2 (miRBase hsa_mature) → UMI deduplication → miRNA counts
```

Supporting modules: `main_downstream.nf` · `main_spikein.nf` · `main_smallrna.nf`

---

## Current Outputs

| Layer | Output | Size / Notes |
|---|---|---|
| cfWGS | 86 aligned BAM files + indexes; 86 MarkDuplicates reports; Picard insert-size reports; Samtools flagstat/idxstats/stats | — |
| cfMeDIP-seq | 98 aligned BAM files + indexes; 300-bp bin count matrix | — |
| Total cfRNA-seq | `totalrna_gene_counts.txt` (protein-coding + lncRNA); `ercc_counts_matrix.tsv` (92 ERCC × 127); `spikein_qc_summary.tsv` | ~52 MB counts; mean ERCC ≈ 2.66%, mean rRNA ≈ 8.66% |
| Small cfRNA-seq | `mirna_counts_matrix.tsv` (1,402 mature miRNAs × 49 samples) | UMI-deduplicated |

---

## Results Showcase

### Upstream Quality Control

![STAR Summary Statistics](docs/images/upstream/star_summary_table%20.png)

![General Statistics](docs/images/upstream/general_stats_table%20.png)

![Picard Insert Size Distribution](docs/images/upstream/picard_insert_size%20.png)

---

## Planned Downstream

Downstream analysis will involve feature extraction, statistical analysis, and feature selection within each omics layer to identify informative molecular features.

- **cfWGS:** Characterization and selection of CNV features, alongside complementary fragmentomic analyses.
- **cfMeDIP-seq:** Identification and selection of informative DNA methylation features.
- **Total cfRNA-seq:** Differential gene-expression analysis and selection of informative transcriptomic features.
- **Small cfRNA-seq:** Characterization and selection of miRNA features, analyzed as a separate branch.

Selected CNV, methylation, and gene-expression features will then be integrated using ActivePathways to investigate cross-omics relationships and enriched biological pathways, followed by network visualization in Cytoscape.

Small cfRNA/miRNA and complementary fragmentomic analyses will remain separate from the core ActivePathways integration.

---

## Repository Structure

```
cellfree-multiomics-gi-cancer/
├── main.nf                        # Nextflow entry point; orchestrates all four branches
├── nextflow.config                # Cluster execution profile (select with -profile cluster)
├── submit_multiomics.sh           # SLURM launcher (1 node, 1 CPU, 16 GB, 120 h)
├── workflows/
│   ├── main_downstream.nf         # Downstream feature generation (CNV, fragmentomics)
│   ├── main_spikein.nf            # ERCC and rRNA spike-in QC module
│   └── main_smallrna.nf           # Small cfRNA / miRNA processing module
├── submit_smallrna.sh             # SLURM launcher for small RNA branch
├── submit_smallrna_array.sh       # Array SLURM launcher for small RNA
├── submit_spikein.sh              # SLURM launcher for spike-in QC module
├── envs/
│   ├── cfwgs.yml                  # Conda environment: cfWGS tools
│   ├── cfmedip.yml                # Conda environment: cfMeDIP-seq tools
│   ├── totalrna.yml               # Conda environment: total cfRNA tools
│   └── smallrna.yml               # Conda environment: small cfRNA / miRNA tools
├── sample_sheets/
│   ├── cfwgs_samplesheet.csv      # SRA run table: cfWGS (86 samples)
│   ├── cfmedip_samplesheet.csv    # SRA run table: cfMeDIP-seq (98 samples)
│   ├── totalrna_samplesheet.csv   # SRA run table: total cfRNA (127 samples)
│   └── smallrna_samplesheet.csv   # SRA run table: small cfRNA (49 samples)
├── scripts/
│   └── umi_dedup.py               # UMI-based deduplication for miRNA counting
├── docs/
│   └── workflow_overview.md       # Detailed methods and parameter notes
└── .gitignore
```

---

## Usage

### Prerequisites

- Nextflow ≥ 22.10
- Conda / Mamba
- SLURM workload manager
- ~2 TB scratch space for all four branches

### Quick start

```bash
# Clone
git clone https://github.com/HABIOTECH/cellfree-multiomics-gi-cancer.git
cd cellfree-multiomics-gi-cancer

# Submit all four branches
sbatch submit_multiomics.sh

# Resume after interruption
nextflow run main.nf -c nextflow.config -profile cluster -resume \
  --cfwgs_samplesheet    sample_sheets/cfwgs_samplesheet.csv \
  --cfmedip_samplesheet  sample_sheets/cfmedip_samplesheet.csv \
  --totalrna_samplesheet sample_sheets/totalrna_samplesheet.csv \
  --smallrna_samplesheet sample_sheets/smallrna_samplesheet.csv
```

---

## Tools and References

| Tool | Version | Use |
|---|---|---|
| Nextflow | DSL2 | Workflow orchestration |
| fasterq-dump (SRA-tools) | — | FASTQ generation from SRA |
| FastQC | — | Per-sample QC |
| Trim Galore | — | Adapter trimming (cfWGS, small cfRNA) |
| fastp | — | Adapter trimming (cfMeDIP-seq) |
| Cutadapt | — | Adapter trimming (total cfRNA) |
| Bowtie2 | — | Alignment: cfWGS, cfMeDIP-seq, ERCC, rRNA, miRNA |
| STAR | — | Splice-aware alignment: total cfRNA (GRCh38) |
| Samtools | — | BAM sorting, indexing, QC |
| GATK MarkDuplicates | — | Duplicate marking |
| Picard | — | Insert-size distribution |
| featureCounts | — | Gene quantification (total cfRNA) |
| miRBase (hsa_mature) | — | miRNA reference (small cfRNA) |
| GRCh38 | 2024-A / GENCODE v44 | Human reference genome |

---

## Human-in-the-Loop AI-Assisted Development

This pipeline was developed with AI assistance (Claude Code) for literature review, Nextflow script drafting, SLURM log interpretation, and workflow debugging. All scientific decisions, HPC execution, output verification, and QC interpretation were performed by the researcher. AI did not run the pipeline or validate results independently.

---

## Citation

If you use this pipeline or its outputs, please cite the original study:

> Tao Y et al. Cell-free multi-omics analysis reveals potential biomarkers in gastrointestinal cancer patients' blood. *Cell Reports Medicine* 2023;4:101281.

---

## Team

**Noura Khaled** — Pipeline development and presentation  
**HOMiCS Bioinformatics Hub** — Horizon Science Communication LLC  
PI: **Hesham M. Abdullah, Ph.D.**

---

## Reference

> Tao Y, Xing S, Zuo S, et al. **Cell-free multi-omics analysis reveals potential biomarkers in gastrointestinal cancer patients' blood.** *Cell Reports Medicine* 2023;4:101281. DOI: [10.1016/j.xcrm.2023.101281](https://doi.org/10.1016/j.xcrm.2023.101281)

## License

MIT License — see [LICENSE](LICENSE) for details.
