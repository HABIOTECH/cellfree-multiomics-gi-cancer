# Cell-Free Multi-Omics Analysis in Gastrointestinal Cancer

**HOMiCS Bioinformatics Hub · Horizon Science Communication LLC**

A production-grade Nextflow DSL2 pipeline reproducing the computational analyses from [Tao et al. 2023 (*Cell Reports Medicine* 4:101281)](https://doi.org/10.1016/j.xcrm.2023.101281) across four public cell-free omics layers in colorectal cancer (CRC), gastric cancer (STAD), and healthy donors.

---

## Reference

> Tao Y, Xing S, Zuo S, et al. **Cell-free multi-omics analysis reveals potential biomarkers in gastrointestinal cancer patients' blood.** *Cell Reports Medicine* 2023;4:101281. DOI: [10.1016/j.xcrm.2023.101281](https://doi.org/10.1016/j.xcrm.2023.101281)

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

## Planned Downstream

Multi-omics integration across cfWGS (CNV features), cfMeDIP-seq (differential methylation), and total cfRNA-seq (differential expression) using **ActivePathways** followed by network visualization in **Cytoscape**. Small cfRNA/miRNA is analyzed separately and is not part of the ActivePathways integration.

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

## License

MIT License — see [LICENSE](LICENSE) for details.
