#!/bin/bash
# =============================================================================
#  spikein_presence_test.sh
#  Cheap go / no-go test: are ERCC and lambda actually present in these
#  libraries?  Subsamples 1 M reads from one sample per layer and maps them
#  against the mini spike-in reference.  Takes a few minutes, no SLURM needed.
#
#  Usage:
#      cd /path/to/cellfree-multiomics-gi-cancer
#      bash spikein_presence_test.sh SRR16575243 SRR16574807
#                                    ^totalRNA    ^cfMeDIP
#
#  If no arguments are given, the first sample of each layer is picked
#  automatically from sample_sheet_CRC_cfOmics.csv.
#
#  INTERPRETATION
#      ERCC   in totalRNA : genuine spike-in is typically 0.5 - 5 %
#                           < 0.01 %  ->  not spiked, drop the ERCC branch
#      lambda in cfMeDIP  : published cfMeDIP protocols target 0.1 - 1 %
#                           < 0.01 %  ->  not spiked, drop the lambda branch
#      rRNA   in totalRNA : expect anywhere from 2 % to 70 %; this is an
#                           UPPER bound (no human genome competing for reads)
# =============================================================================
set -euo pipefail

# ── User-configurable paths ────────────────────────────────────────────────────
PROJ="${PROJECT_DIR:-$(pwd)}"
# ──────────────────────────────────────────────────────────────────────────────

SPK="${PROJ}/references/spikein"
IDX="${SPK}/spikein"
SHEET="${PROJ}/sample_sheet_CRC_cfOmics.csv"
OUT="${PROJ}/results/spikein_presence_test"
NREADS=1000000          # reads to subsample (not lines)

BOWTIE2_MOD="bowtie2/2.2.7"
SAMTOOLS_MOD="samtools/1.21"

cd "${PROJ}"
mkdir -p "${OUT}"

if [[ ! -s "${IDX}.1.bt2" ]]; then
    echo "ERROR: bowtie2 index not found at ${IDX}. Run build_spikein_refs.sh first."
    exit 1
fi

# ------------------------------------------------------------- pick samples
pick_layer () {   # $1 = layer name
    awk -F, -v L="$1" 'NR>1 && $2==L {print $1; exit}' "${SHEET}"
}
RNA_ID="${1:-$(pick_layer totalRNA)}"
MED_ID="${2:-$(pick_layer cfMeDIP)}"

# --------------------------------------------- locate the trimmed R1 FASTQ
# trimgalore : SRR*_1_val_1.fq.gz      (cfDNA + smallRNA)
# fastp      : SRR*_trimmed_1.fastq.gz (cfMeDIP)
# cutadapt   : SRR*_trimmed_1.fastq.gz (totalRNA)
find_r1 () {      # $1 = SRR id
    local id="$1" f
    for f in "results/trimmed/cutadapt/${id}_trimmed_1.fastq.gz" \
             "results/trimmed/fastp/${id}_trimmed_1.fastq.gz" \
             "results/trimmed/trimgalore/${id}_1_val_1.fq.gz"; do
        [[ -s "$f" ]] && { echo "$f"; return 0; }
    done
    # fall back to anything inside work/
    f=$(find work -name "${id}*_val_1.fq.gz" -o -name "${id}_trimmed_1.fastq.gz" 2>/dev/null | head -1)
    [[ -n "$f" ]] && { echo "$f"; return 0; }
    return 1
}

module load ${BOWTIE2_MOD}
module load ${SAMTOOLS_MOD}

run_one () {      # $1 = SRR id, $2 = label
    local id="$1" label="$2" fq
    echo
    echo "=============================================================="
    echo " ${label}  :  ${id}"
    echo "=============================================================="
    if ! fq=$(find_r1 "${id}"); then
        echo "  ERROR: no trimmed R1 FASTQ found for ${id}"
        return
    fi
    echo "  FASTQ : ${fq}"

    # subsample without seqtk: take the first NREADS records
    zcat "${fq}" 2>/dev/null | head -n $((NREADS * 4)) > "${OUT}/${id}.sub.fq" || true
    local n
    n=$(( $(wc -l < "${OUT}/${id}.sub.fq") / 4 ))
    echo "  reads : ${n}"

    bowtie2 -x "${IDX}" -U "${OUT}/${id}.sub.fq" -p 4 --no-unal \
        2> "${OUT}/${id}.bowtie2.log" \
      | samtools sort -o "${OUT}/${id}.spikein.bam" -
    samtools index "${OUT}/${id}.spikein.bam"
    rm -f "${OUT}/${id}.sub.fq"

    # aggregate idxstats by class
    samtools idxstats "${OUT}/${id}.spikein.bam" \
      | awk -v tot="${n}" -v map="${SPK}/contig_map.tsv" '
            BEGIN{ OFS="\t"; while ((getline line < map) > 0){ split(line,a,"\t"); cls[a[1]]=a[2] } }
            $1!="*" { c = (cls[$1] ? cls[$1] : "other"); s[c] += $3 }
            END {
                printf "  %-14s %12s %10s\n", "CLASS", "READS", "PCT"
                for (k in s) printf "  %-14s %12d %9.4f%%\n", k, s[k], 100*s[k]/tot
            }' | sort

    echo "  --- bowtie2 overall ---"
    grep "overall alignment rate" "${OUT}/${id}.bowtie2.log" | sed 's/^/  /'
}

run_one "${RNA_ID}" "totalRNA (ERCC + rRNA expected)"
run_one "${MED_ID}" "cfMeDIP  (lambda expected)"

echo
echo "=============================================================="
echo " Results written to ${OUT}"
echo "=============================================================="
