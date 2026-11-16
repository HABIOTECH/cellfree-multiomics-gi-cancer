#!/bin/bash
# =============================================================================
#  build_spikein_refs.sh
#  Builds a lightweight spike-in / contaminant reference panel (Option B).
#
#  Run from the project root:
#      cd /path/to/cellfree-multiomics-gi-cancer
#      bash build_spikein_refs.sh
#
#  PREREQUISITE — upload these via FileZilla (or scp) into references/spikein/raw/ :
#      ERCC92.fa            (Thermo Fisher ERCC control analysis file)
#      cms_095046.txt       (ERCC nominal concentration table)
#      lambda.fa            (NCBI NC_001416.1, 48502 bp)
#      rDNA_U13369.fa       (NCBI U13369.1, human rDNA repeat unit, 42999 bp)
#      UniVec_Core          (NCBI ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/)
#      [optional] medip_spikein.fa   (paper's synthetic me/unme amplicons)
#
#  Produces:
#      references/spikein/spikein.fa          concatenated mini-reference
#      references/spikein/spikein.*.bt2       bowtie2 index
#      references/spikein/contig_map.tsv      contig -> class lookup
#      references/spikein/univec/             blast db for VecScreen
# =============================================================================
set -euo pipefail

# ── User-configurable paths ────────────────────────────────────────────────────
PROJ="${PROJECT_DIR:-$(pwd)}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
# ──────────────────────────────────────────────────────────────────────────────

SPK="${PROJ}/references/spikein"
RAW="${SPK}/raw"

BOWTIE2_MOD="bowtie2/2.2.7"
SAMTOOLS_MOD="samtools/1.21"

cd "${PROJ}"
mkdir -p "${RAW}" "${SPK}/univec"

echo "=============================================================="
echo " STEP 1 : check uploaded reference files"
echo "=============================================================="
missing=0
for f in ERCC92.fa lambda.fa rDNA_U13369.fa; do
    if [[ -s "${RAW}/${f}" ]]; then
        n=$(grep -c '^>' "${RAW}/${f}")
        printf "  OK       %-20s %6d sequence(s)\n" "$f" "$n"
    else
        printf "  MISSING  %-20s  <-- upload to %s\n" "$f" "${RAW}"
        missing=1
    fi
done
for f in cms_095046.txt UniVec_Core medip_spikein.fa; do
    if [[ -s "${RAW}/${f}" ]]; then
        printf "  OK       %-20s (optional)\n" "$f"
    else
        printf "  absent   %-20s (optional)\n" "$f"
    fi
done
if [[ $missing -eq 1 ]]; then
    echo
    echo "ERROR: required reference(s) missing. Upload them and re-run."
    exit 1
fi

echo
echo "=============================================================="
echo " STEP 2 : normalise contig names and concatenate"
echo "=============================================================="
awk '/^>/{print $1; next}{print}' "${RAW}/ERCC92.fa" > "${SPK}/.ercc.fa"
awk 'BEGIN{n=0} /^>/{n++; print (n==1 ? ">lambda" : ">lambda_"n); next}{print}' \
    "${RAW}/lambda.fa" > "${SPK}/.lambda.fa"
awk 'BEGIN{n=0} /^>/{n++; print (n==1 ? ">rDNA_U13369" : ">rDNA_"n); next}{print}' \
    "${RAW}/rDNA_U13369.fa" > "${SPK}/.rdna.fa"

cat "${SPK}/.ercc.fa" "${SPK}/.lambda.fa" "${SPK}/.rdna.fa" > "${SPK}/spikein.fa"

if [[ -s "${RAW}/medip_spikein.fa" ]]; then
    awk '/^>/{print $1; next}{print}' "${RAW}/medip_spikein.fa" >> "${SPK}/spikein.fa"
    echo "  medip_spikein.fa appended"
fi
rm -f "${SPK}/.ercc.fa" "${SPK}/.lambda.fa" "${SPK}/.rdna.fa"

echo "  spikein.fa : $(grep -c '^>' "${SPK}/spikein.fa") contigs, \
$(grep -v '^>' "${SPK}/spikein.fa" | tr -d '\n' | wc -c) bp"

echo
echo "=============================================================="
echo " STEP 3 : contig -> class lookup table"
echo "=============================================================="
grep '^>' "${SPK}/spikein.fa" | sed 's/^>//' | awk 'BEGIN{OFS="\t"}
{
    if ($1 ~ /^ERCC-/)        cls="ERCC"
    else if ($1 ~ /^lambda/)  cls="lambda"
    else if ($1 ~ /^rDNA/)    cls="rRNA"
    else                      cls="medip_spike"
    print $1, cls
}' > "${SPK}/contig_map.tsv"
cut -f2 "${SPK}/contig_map.tsv" | sort | uniq -c | sed 's/^/  /'

echo
echo "=============================================================="
echo " STEP 4 : build bowtie2 index"
echo "=============================================================="
module load ${BOWTIE2_MOD}
module load ${SAMTOOLS_MOD}
bowtie2-build --threads 4 "${SPK}/spikein.fa" "${SPK}/spikein" > "${SPK}/bowtie2-build.log" 2>&1
samtools faidx "${SPK}/spikein.fa"
ls -lh "${SPK}"/spikein*.bt2 | sed 's/^/  /'

echo
echo "=============================================================="
echo " STEP 5 : UniVec blast database (optional)"
echo "=============================================================="
if [[ -s "${RAW}/UniVec_Core" ]]; then
    source "${CONDA_SH}"
    conda activate blast_env
    cp "${RAW}/UniVec_Core" "${SPK}/univec/UniVec_Core.fa"
    makeblastdb -in "${SPK}/univec/UniVec_Core.fa" -dbtype nucl \
                -out "${SPK}/univec/UniVec_Core" > "${SPK}/univec/makeblastdb.log" 2>&1
    echo "  UniVec_Core blast db built"
    conda deactivate
else
    echo "  skipped (UniVec_Core not uploaded)"
fi

echo
echo "=============================================================="
echo " DONE"
echo "=============================================================="
echo "  Index prefix : ${SPK}/spikein"
echo "  Next step    : bash spikein_presence_test.sh <totalRNA_SRR> <cfMeDIP_SRR>"
