nextflow.enable.dsl=2

workflow {
    dna_ids       = Channel.fromPath(params.dna_csv).splitCsv(header:true).map{ row -> tuple(row.Run, 'cfDNA') }
    medip_ids     = Channel.fromPath(params.medip_csv).splitCsv(header:true).map{ row -> tuple(row.Run, 'cfMeDIP') }
    total_rna_ids = Channel.fromPath(params.total_rna_csv).splitCsv(header:true).map{ row -> tuple(row.Run, 'totalRNA') }
    small_rna_ids = Channel.fromPath(params.small_rna_csv).splitCsv(header:true).map{ row -> tuple(row.Run, 'smallRNA') }

    all_ids = dna_ids.concat(medip_ids, total_rna_ids, small_rna_ids)

    // 1. Build DNA Index if it doesn't exist
    INDEX_DNA(params.genome_fasta)

    // 2. Download and QC
    DOWNLOAD_SRA(all_ids)
    FASTQC_RAW(DOWNLOAD_SRA.out.fastqs)

    // 3. Trimming — tool chosen per omics layer
    DOWNLOAD_SRA.out.fastqs.branch {
        trimgalore: it[1] == 'cfDNA' || it[1] == 'smallRNA'
        fastp:      it[1] == 'cfMeDIP'
        cutadapt:   it[1] == 'totalRNA'
    }.set { branched_raw }

    TRIM_GALORE(branched_raw.trimgalore)
    TRIM_FASTP(branched_raw.fastp)
    TRIM_CUTADAPT(branched_raw.cutadapt)

    trimmed_all = TRIM_GALORE.out.trimmed.concat(TRIM_FASTP.out.trimmed, TRIM_CUTADAPT.out.trimmed)
    FASTQC_TRIM(trimmed_all)

    // 4. Alignment — BWA/Bowtie2 for DNA layers, STAR for RNA layers
    trimmed_all.branch {
        dna: it[1] == 'cfDNA' || it[1] == 'cfMeDIP'
        rna: it[1] == 'totalRNA' || it[1] == 'smallRNA'
    }.set { branched_trimmed }

    ALIGN_DNA(branched_trimmed.dna, INDEX_DNA.out.index.collect())
    ALIGN_RNA(branched_trimmed.rna)
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESSES
// ─────────────────────────────────────────────────────────────────────────────

process INDEX_DNA {
    tag "Bowtie2_Index"
    publishDir "${params.outdir}/references/bowtie2_index", mode: 'copy'
    input:  path fasta
    output: path "genome*", emit: index
    script:
    """
    module load ${params.modules.bowtie2}
    bowtie2-build --threads ${task.cpus} ${fasta} genome
    """
}

process DOWNLOAD_SRA {
    tag "${id}"
    input:  tuple val(id), val(type)
    output: tuple val(id), val(type), path("${id}_{1,2}.fastq.gz"), emit: fastqs
    script:
    """
    module load ${params.modules.sra}
    fasterq-dump ${id} --split-files
    gzip ${id}_1.fastq
    gzip ${id}_2.fastq
    """
}

process FASTQC_RAW {
    tag "${id}"
    input:  tuple val(id), val(type), path(reads)
    output: path "*.zip", emit: zips
    script:
    """
    module load ${params.modules.fastqc}
    fastqc -t ${task.cpus} ${reads}
    """
}

process TRIM_GALORE {
    tag "${id}"
    publishDir "${params.outdir}/trimmed/trimgalore", mode: 'copy'
    input:  tuple val(id), val(type), path(reads)
    output: tuple val(id), val(type), path("${id}*_val_{1,2}.fq.gz"), emit: trimmed
    script:
    """
    # Activate trim-galore conda environment
    # Set CONDA_SH to your conda initialisation script, e.g.:
    #   export CONDA_SH=\${HOME}/miniconda3/etc/profile.d/conda.sh
    source "\${CONDA_SH:-\${HOME}/miniconda3/etc/profile.d/conda.sh}"
    conda activate trimgalore
    trim_galore --cores 4 --paired --gzip ${reads[0]} ${reads[1]}
    """
}

process TRIM_FASTP {
    tag "${id}"
    publishDir "${params.outdir}/trimmed/fastp", mode: 'copy'
    input:  tuple val(id), val(type), path(reads)
    output: tuple val(id), val(type), path("${id}_trimmed_{1,2}.fastq.gz"), emit: trimmed
    script:
    """
    module load ${params.modules.fastp}
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${id}_trimmed_1.fastq.gz -O ${id}_trimmed_2.fastq.gz
    """
}

process TRIM_CUTADAPT {
    tag "${id}"
    publishDir "${params.outdir}/trimmed/cutadapt", mode: 'copy'
    input:  tuple val(id), val(type), path(reads)
    output: tuple val(id), val(type), path("${id}_trimmed_{1,2}.fastq.gz"), emit: trimmed
    script:
    """
    module load ${params.modules.cutadapt}
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -m 20 \
        -o ${id}_trimmed_1.fastq.gz -p ${id}_trimmed_2.fastq.gz \
        ${reads[0]} ${reads[1]}
    """
}

process FASTQC_TRIM {
    tag "${id}"
    input:  tuple val(id), val(type), path(reads)
    output: path "*.zip", emit: zips
    script:
    """
    module load ${params.modules.fastqc}
    fastqc -t ${task.cpus} ${reads}
    """
}

process ALIGN_DNA {
    tag "${id}"
    publishDir "${params.outdir}/alignment/dna", mode: 'copy'
    input:
        tuple val(id), val(type), path(reads)
        path index
    output: tuple val(id), val(type), path("${id}.bam"), emit: bam
    script:
    """
    module load ${params.modules.bowtie2}
    module load ${params.modules.samtools}
    bowtie2 -x genome -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} | \
    samtools view -bS - | samtools sort -o ${id}.bam
    """
}

process ALIGN_RNA {
    tag "${id}"
    publishDir "${params.outdir}/alignment/rna", mode: 'copy'
    input:  tuple val(id), val(type), path(reads)
    output: tuple val(id), val(type), path("${id}_Aligned.sortedByCoord.out.bam"), emit: bam
    script:
    """
    module load ${params.modules.star}
    STAR --genomeDir ${params.star_index} \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${id}_ \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 30000000000 \
         --runThreadN ${task.cpus}
    """
}
