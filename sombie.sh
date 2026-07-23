#!/bin/bash
#
#SBATCH -J sombie
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 3-00:00:00

# --- Usage function ---
usage() {
  cat <<EOF
This script creates a filtered somatic VCF file from matched tumor/normal or tumor-only BAM files using GATK Mutect2.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] sombie.sh -I <tumor.bam> [-N <normal.bam>] -O <out_dir> -R <ref.fa> [-L <intervals.bed>] ...

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
By default the script will use 32 threads to utilize per chromosome parallelization and speed up variant calling.
Do not change the number of threads unless using intervals (-L) or in --singlethread mode.
The best way to name log files is sample name followed by job name and ID (-o  .../logs/SAMPLE_%x_%A.log).

Mandatory flags:
  -I <file>      Path to the tumor BAM file. You can use the output of bampire.sh.
  -O <dir>       Path to the desired output directory.
  -R <file>      Path to the reference genome FASTA file to which the BAM files are aligned.

Optional flags:
  -N <file>                  Path to the matched normal BAM file. When this is added the calling switches to tumor-normal.
  -S <name>                  Tumor sample name, used for naming output files. If not provided, it is read from the tumor BAM's read group (SM tag).
  -NS <name>                 Normal sample name, passed to Mutect2 as the -normal argument. If not provided, it is read from the normal BAM's read group (SM tag).
  -L <file>                  Path to a BED file with target intervals for exome or targeted sequencing.
                             This disables multi-threading so make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.
  --singlethread             Disable per chromosome parallelization for genome sequencing and perform variant calling on a single cpu thread.
                             This will be significantly slower but can be used when running a large number of samples in parallel (>50).
                             Make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.
  --custom-pon <file>        Path to a custom Panel of Normals VCF.
  --custom-germline <file>   Path to a custom germline resource VCF.
  --custom-common <file>     Path to a custom small biallelic common-SNP VCF for contamination estimation.
  --genotype-pon-sites       Genotype and filter out sites that match the panel of normals instead of dropping them.
  --genotype-germline-sites  Genotype and filter out sites that look germline per the germline resource instead of dropping them.
  --normal-lod <float>       Override Mutect2's --normal-lod threshold (default 2.2).
                             Higher values make it harder to reject a candidate based on normal-sample evidence (more permissive toward calling somatic).
  --disable-mate-filter      Disable MateOnSameContigOrNoMappedMateReadFilter, retaining read pairs whose mate maps to a different contig.
                             Useful for capturing evidence near translocation breakpoints; adds noise.
  --min-depth <int>          Filter variants whose tumor total depth is below this threshold.
  --min-alt-reads <int>      Filter variants whose tumor ALT allele count is below this threshold.
  --min-vaf                  Filter sites below a certain VAF in the tumor (value range from 0 to 1, eg. 0.10 for 10%)
  --blacklist-filter         Filter variants overlapping the Encode blacklist region file.
  --custom-blacklist <file>  Path to a custom BED file of blacklist/low-complexity regions.


EOF
  exit 1
}

# --- Default resources ---
DEFAULT_PON="/lustre/imgge/lab01/refs/db/hg38/gatk_1000g_pon.hg38.vcf.gz"
DEFAULT_GERMLINE="/lustre/imgge/lab01/refs/db/hg38/gatk_af-only-gnomad.hg38.vcf.gz"
DEFAULT_COMMON="/lustre/imgge/lab01/refs/db/hg38/gatk_af-only-gnomad-common-biallelic.vcf.gz"
# DEFAULT_BWA_INDEX_IMAGE="/lustre/imgge/lab01/refs/hg38_gatk/hg38.fasta.img"
DEFAULT_BLACKLIST="/lustre/imgge/lab01/refs/db/hg38/hg38_blacklist.v2.bed.gz"

# --- Initial check ---
if [ "$#" -eq 0 ]; then
  usage
fi

set -euo pipefail

# --- Argument Parsing ---
TUMOR_BAM=""
OUTPUT_DIR=""
REF=""
NORMAL_BAM=""
TUMOR_SAMPLE=""
NORMAL_SAMPLE=""
INTERVAL_FILE=""
INTERVALS=""
PON="${DEFAULT_PON}"
GERMLINE_RESOURCE="${DEFAULT_GERMLINE}"
CONTAM_RESOURCE="${DEFAULT_COMMON}"
# BWA_INDEX_IMAGE="${DEFAULT_BWA_INDEX_IMAGE}"
# RUN_ALIGNMENT_ARTIFACTS=false
BLACKLIST_FILE="${DEFAULT_BLACKLIST}"
GENOTYPE_PON=false
GENOTYPE_GERMLINE=false
NORMAL_LOD=""
DISABLE_MATE_FILTER=false
MIN_DEPTH=""
MIN_ALT_READS=""
MIN_VAF=""
RUN_BLACKLIST=false
SINGLETHREAD_MODE=false

# Manual loop to process options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -I) TUMOR_BAM="$2"; shift 2 ;;
    -N) NORMAL_BAM="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    --custom-germline) GERMLINE_RESOURCE="$2"; shift 2 ;;
    --custom-pon) PON="$2"; shift 2 ;;
    --custom-common) CONTAM_RESOURCE="$2"; shift 2 ;;
    # --custom-bwa-index-image) BWA_INDEX_IMAGE="$2"; shift 2 ;;
    # --filter-alignment) RUN_ALIGNMENT_ARTIFACTS=true; shift ;;
    -S) TUMOR_SAMPLE="$2"; shift 2 ;;
    -NS) NORMAL_SAMPLE="$2"; shift 2 ;;
    -L)
      INTERVAL_FILE="$2"
      if [[ "${INTERVAL_FILE}" != *.bed ]]; then echo "Error: -L file must be a .bed file." >&2; usage; fi
      if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: ${INTERVAL_FILE}" >&2; usage; fi
      INTERVALS="-L ${INTERVAL_FILE} -ip 20"
      shift 2
      ;;
    --genotype-pon-sites) GENOTYPE_PON=true; shift ;;
    --genotype-germline-sites) GENOTYPE_GERMLINE=true; shift ;;
    --normal-lod) NORMAL_LOD="$2"; shift 2 ;;
    --disable-mate-filter) DISABLE_MATE_FILTER=true; shift ;;
    --singlethread) SINGLETHREAD_MODE=true; shift ;;
    --min-depth) MIN_DEPTH="$2"; shift 2 ;;
    --min-alt-reads) MIN_ALT_READS="$2"; shift 2 ;;
    --min-vaf) MIN_VAF="$2"; shift 2 ;;
    --blacklist-filter) RUN_BLACKLIST=true; shift ;;
    --custom-blacklist) BLACKLIST_FILE="$2"; shift 2 ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$TUMOR_BAM" ]]; then echo "Error: -I <tumor.bam> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi
if [[ ! -f "$GERMLINE_RESOURCE" ]]; then
  echo "Error: Germline resource not found at ${GERMLINE_RESOURCE}. Provide --custom-germline <file>." >&2
  exit 1
fi
if [[ ! -f "$PON" ]]; then
  echo "Error: Panel of Normals not found at ${PON}. Provide --custom-pon <file>." >&2
  exit 1
fi
if [[ ! -f "$CONTAM_RESOURCE" ]]; then
  echo "Error: Common-SNP contamination resource not found at ${CONTAM_RESOURCE}. Provide --custom-common <file>." >&2
  exit 1
fi
# if [ "$RUN_ALIGNMENT_ARTIFACTS" = true ] && [[ ! -f "$BWA_INDEX_IMAGE" ]]; then
#   echo "Error: BWA index image not found at ${BWA_INDEX_IMAGE}. Provide --custom-bwa-index-image <file>." >&2
#   exit 1
# fi
if [ "$RUN_BLACKLIST" = true ] && [[ ! -f "$BLACKLIST_FILE" ]]; then
  echo "Error: Blacklist file not found at ${BLACKLIST_FILE}. Provide --custom-blacklist <file>." >&2
  exit 1
fi
if [[ -n "$MIN_DEPTH" ]] && ! [[ "$MIN_DEPTH" =~ ^[0-9]+$ ]]; then
  echo "Error: --min-depth must be a non-negative integer (e.g. 10)." >&2
  exit 1
fi
if [[ -n "$MIN_ALT_READS" ]] && ! [[ "$MIN_ALT_READS" =~ ^[0-9]+$ ]]; then
  echo "Error: --min-alt-reads must be a non-negative integer (e.g. 3)." >&2
  exit 1
fi
if [[ -n "$MIN_VAF" ]] && ! [[ "$MIN_VAF" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]]; then
  echo "Error: --min-vaf must be a number between 0 and 1 (e.g. 0.01 for 1%)." >&2
  exit 1
fi


# --- Start script ---
# Create output directories
mkdir -p "${OUTPUT_DIR}/vcfs/metrics"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

# Derive sample names from BAM read groups if not provided
if [[ -z "$TUMOR_SAMPLE" ]]; then
  echo "INFO: Tumor sample name not provided with -S. Deriving from tumor BAM read group (SM tag)."
  TUMOR_SAMPLE=$(samtools view -H "${TUMOR_BAM}" | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/' | head -1)
  if [[ -z "$TUMOR_SAMPLE" ]]; then echo "Error: Could not derive tumor sample name from BAM read group." >&2; exit 1; fi
fi

if [[ -n "$NORMAL_BAM" ]]; then
  TUMOR_ONLY=false
  if [[ -z "$NORMAL_SAMPLE" ]]; then
    echo "INFO: Normal sample name not provided with -NS. Deriving from normal BAM read group (SM tag)."
    NORMAL_SAMPLE=$(samtools view -H "${NORMAL_BAM}" | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/' | head -1)
    if [[ -z "$NORMAL_SAMPLE" ]]; then echo "Error: Could not derive normal sample name from BAM read group." >&2; exit 1; fi
  fi
  echo "INFO: Tumor sample: ${TUMOR_SAMPLE} | Normal sample: ${NORMAL_SAMPLE}"
  PREFIX="${TUMOR_SAMPLE}-tumor-normal"
else
  TUMOR_ONLY=true
  echo "INFO: No -N <normal.bam> provided. Running in TUMOR-ONLY mode for sample: ${TUMOR_SAMPLE}"
  if [ "$GENOTYPE_GERMLINE" = false ]; then
    echo "WARNING: Running tumor-only without --genotype-germline-sites. Without a matched normal, germline variants" >&2
    echo "         are rejected using only the germline resource/PoN. Consider adding --genotype-germline-sites so" >&2
    echo "         germline-flagged sites remain visible (tagged, not dropped) for manual review." >&2
  fi
  PREFIX="${TUMOR_SAMPLE}-tumor-only"
fi

# Build optional Mutect2 argument string
NORMAL_ARGS=()
if [ "$TUMOR_ONLY" = false ]; then
  NORMAL_ARGS=("-I" "${NORMAL_BAM}" "-normal" "${NORMAL_SAMPLE}")
fi

MUTECT_EXTRA_ARGS=()
if [ "$GENOTYPE_PON" = true ]; then
  MUTECT_EXTRA_ARGS+=("--genotype-pon-sites")
fi
if [ "$GENOTYPE_GERMLINE" = true ]; then
  MUTECT_EXTRA_ARGS+=("--genotype-germline-sites")
fi
if [[ -n "$NORMAL_LOD" ]]; then
  MUTECT_EXTRA_ARGS+=("--normal-lod" "${NORMAL_LOD}")
fi
if [ "$DISABLE_MATE_FILTER" = true ]; then
  MUTECT_EXTRA_ARGS+=("--disable-read-filter" "MateOnSameContigOrNoMappedMateReadFilter")
fi

# --- Somatic Variant Calling ---
# Run Mutect2
if [[ "$SINGLETHREAD_MODE" = true || -n "$INTERVAL_FILE" ]]; then
  echo "INFO: Running Mutect2 in single-thread mode..."
  gatk Mutect2 \
    -R "${REF}" \
    ${INTERVALS} \
    -I "${TUMOR_BAM}" \
    "${NORMAL_ARGS[@]}" \
    --germline-resource "${GERMLINE_RESOURCE}" \
    --panel-of-normals "${PON}" \
    --f1r2-tar-gz "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz" \
    "${MUTECT_EXTRA_ARGS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${PREFIX}_raw.vcf.gz"
  echo "INFO: Finished Mutect2 in single-thread mode!"
 
  F1R2_ARGS=("-I" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz")
  mv "${OUTPUT_DIR}/vcfs/${PREFIX}_raw.vcf.gz.stats" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_mutect.stats"
else
  CHRS=(chr{1..22} chrX chrY chrM)
  CHR_VCFS=()
  STATS_ARGS=()
  F1R2_ARGS=()
  echo "INFO: Running Mutect2 per chromosome..."
  for CHR in "${CHRS[@]}"; do
    CHR_VCFS+=("-I" "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}.vcf.gz")
    STATS_ARGS+=("-stats" "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}.vcf.gz.stats")
    F1R2_ARGS+=("-I" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_${CHR}_f1r2.tar.gz")
    (
      gatk Mutect2 \
        -R "${REF}" \
        -L "${CHR}" \
        -I "${TUMOR_BAM}" \
        "${NORMAL_ARGS[@]}" \
        --germline-resource "${GERMLINE_RESOURCE}" \
        --panel-of-normals "${PON}" \
        --f1r2-tar-gz "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_${CHR}_f1r2.tar.gz" \
        "${MUTECT_EXTRA_ARGS[@]}" \
        -O "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}.vcf.gz"
      echo "INFO: Finished Mutect2 for ${CHR}!"
    ) &
  done
  wait
  echo "INFO: All chromosome jobs finished!"
 
  echo "INFO: Merging per-chromosome VCFs..."
  gatk MergeVcfs \
    "${CHR_VCFS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${PREFIX}_raw.vcf.gz"
  echo "INFO: Successfully merged VCF!"
 
  echo "INFO: Merging per-chromosome stats files..."
  gatk MergeMutectStats \
    "${STATS_ARGS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_mutect.stats"
  echo "INFO: Successfully merged stats!"
fi

echo "INFO: Learning read orientation model from F1R2 counts..."
gatk LearnReadOrientationModel \
  "${F1R2_ARGS[@]}" \
  -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_read-orientation-model.tar.gz"
echo "INFO: Finished learning read orientation model!"

# --- Contamination Estimation ---
FILTER_EXTRA_ARGS=()
echo "INFO: Running GetPileupSummaries for tumor..."
gatk GetPileupSummaries \
  -I "${TUMOR_BAM}" \
  -V "${CONTAM_RESOURCE}" \
  -L "${CONTAM_RESOURCE}" \
  ${INTERVALS} \
  -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table"
 
if [ "$TUMOR_ONLY" = false ]; then
  echo "INFO: Running GetPileupSummaries for normal..."
  gatk GetPileupSummaries \
    -I "${NORMAL_BAM}" \
    -V "${CONTAM_RESOURCE}" \
    -L "${CONTAM_RESOURCE}" \
    ${INTERVALS} \
    -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table"
 
  echo "INFO: Calculating contamination (matched normal)..."
  gatk CalculateContamination \
    -I "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table" \
    -matched "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table" \
    -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_contamination.table" \
    --tumor-segmentation "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_segments.table"
else
  echo "INFO: Calculating contamination (tumor-only, no matched normal)..."
  gatk CalculateContamination \
    -I "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table" \
    -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_contamination.table" \
    --tumor-segmentation "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_segments.table"
fi
echo "INFO: Finished contamination estimation!"
 
FILTER_EXTRA_ARGS+=("--contamination-table" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_contamination.table")
FILTER_EXTRA_ARGS+=("--tumor-segmentation" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_segments.table")


# --- Variant Filtering ---
echo "INFO: Filtering VCF..."
gatk FilterMutectCalls \
  -R "${REF}" \
  -V "${OUTPUT_DIR}/vcfs/${PREFIX}_raw.vcf.gz" \
  --ob-priors "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_read-orientation-model.tar.gz" \
  --stats "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_mutect.stats" \
  "${FILTER_EXTRA_ARGS[@]}" \
  -O "${OUTPUT_DIR}/vcfs/${PREFIX}_filter1.vcf.gz"
CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter1.vcf.gz"
echo "INFO: Finished filtering VCF!"

 
# if [ "$RUN_ALIGNMENT_ARTIFACTS" = true ]; then
#   echo "INFO: Running FilterAlignmentArtifacts..."
#   if [[ "$SINGLETHREAD_MODE" = true || -n "$INTERVAL_FILE" ]]; then
#     gatk FilterAlignmentArtifacts \
#       -R "${REF}" \
#       -V "${OUTPUT_DIR}/vcfs/${PREFIX}_filter1.vcf.gz" \
#       -I "${TUMOR_BAM}" \
#       --bwa-mem-index-image "${BWA_INDEX_IMAGE}" \
#       -O "${OUTPUT_DIR}/vcfs/${PREFIX}_filter2.vcf.gz"
#   else
#     ARTIFACT_CHR_VCFS=()
#     for CHR in "${CHRS[@]}"; do
#       ARTIFACT_CHR_VCFS+=("-I" "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}_artifact.vcf.gz")
#       (
#         gatk FilterAlignmentArtifacts \
#           -R "${REF}" \
#           -L "${CHR}" \
#           -V "${OUTPUT_DIR}/vcfs/${PREFIX}_filter1.vcf.gz" \
#           -I "${TUMOR_BAM}" \
#           --bwa-mem-index-image "${BWA_INDEX_IMAGE}" \
#           -O "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}_artifact.vcf.gz"
#         echo "INFO: Finished FilterAlignmentArtifacts for ${CHR}!"
#       ) &
#     done
#     wait
#     echo "INFO: All FilterAlignmentArtifacts chromosome jobs finished!"
#     gatk MergeVcfs \
#       "${ARTIFACT_CHR_VCFS[@]}" \
#       -O "${OUTPUT_DIR}/vcfs/${PREFIX}_filter2.vcf.gz"
#   fi
#   echo "INFO: Finished FilterAlignmentArtifacts!"
#   CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter2.vcf.gz"
# else
#   echo "INFO: Skipping FilterAlignmentArtifacts (--filter-alignment not set)."
#   CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter1.vcf.gz"
# fi

if [ "$TUMOR_ONLY" = false ]; then
  FIRST_SAMPLE=$(bcftools query -l "${CURRENT_VCF}" | head -1)
  if [[ "$FIRST_SAMPLE" = "$NORMAL_SAMPLE" ]]; then
    echo "INFO: Tumor sample is not first (found '${FIRST_SAMPLE}'). Reordering sample columns (tumor first, normal second)..."
    TRUE_TUMOR_SAMPLE=$(bcftools query -l "${CURRENT_VCF}" | head -2 | tail -1)
    bcftools view -s "${TRUE_TUMOR_SAMPLE},${NORMAL_SAMPLE}" \
      "${CURRENT_VCF}" \
      -O z \
      -o "${OUTPUT_DIR}/vcfs/${PREFIX}_reordered.vcf.gz"
    CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_reordered.vcf.gz"
    echo "INFO: Finished reordering sample columns!"
  else
    echo "INFO: Tumor sample is already first in the VCF. Skipping reorder."
  fi
fi

if [[ -n "$MIN_DEPTH" ]]; then
  echo "INFO: Tagging variants with tumor total depth < ${MIN_DEPTH} as LowDepth in FILTER"
  bcftools filter \
    -e "FORMAT/DP[0]<${MIN_DEPTH}" \
    -s "LowDepth" \
    -m+ \
    "${CURRENT_VCF}" \
    -O z \
    -o "${OUTPUT_DIR}/vcfs/${PREFIX}_filter3.vcf.gz"
  tabix "${OUTPUT_DIR}/vcfs/${PREFIX}_filter3.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter3.vcf.gz"
  echo "INFO: Finished LowDepth tagging!"
fi
 
if [[ -n "$MIN_ALT_READS" ]]; then
  echo "INFO: Tagging variants with tumor ALT read count < ${MIN_ALT_READS} as LowAltReads in FILTER"
  bcftools filter \
    -e "FORMAT/AD[0:1]<${MIN_ALT_READS}" \
    -s "LowAltReads" \
    -m+ \
    "${CURRENT_VCF}" \
    -O z \
    -o "${OUTPUT_DIR}/vcfs/${PREFIX}_filter4.vcf.gz"
  tabix "${OUTPUT_DIR}/vcfs/${PREFIX}_filter4.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter4.vcf.gz"
  echo "INFO: Finished LowAltReads tagging!"
fi

if [[ -n "$MIN_VAF" ]]; then
  echo "INFO: Tagging variants with tumor AF < ${MIN_VAF} as LowVAF in FILTER"
  bcftools filter \
    -e "FORMAT/AF[0:0]<${MIN_VAF}" \
    -s "LowVAF" \
    -m+ \
    "${CURRENT_VCF}" \
    -O z \
    -o "${OUTPUT_DIR}/vcfs/${PREFIX}_filter5.vcf.gz"
  tabix "${OUTPUT_DIR}/vcfs/${PREFIX}_filter5.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter5.vcf.gz"
  echo "INFO: Finished LowVAF tagging!"
fi
 
if [ "$RUN_BLACKLIST" = true ]; then
  echo "INFO: Tagging variants overlapping blacklist regions..."
  echo '##FILTER=<ID=Blacklist,Description="Overlaps a known blacklist/low-complexity region (ENCODE blacklist)">' > "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_blacklist_header.txt"
  bcftools annotate \
    -a "${BLACKLIST_FILE}" \
    -c CHROM,FROM,TO \
    --mark-sites "+BLACKLIST_REGION" \
    "${CURRENT_VCF}" \
    -O u | \
  bcftools filter \
    -i "BLACKLIST_REGION=1" \
    -s "Blacklist" \
    -m+ \
    -O z -o "${OUTPUT_DIR}/vcfs/${PREFIX}_filter6.vcf.gz"
  tabix "${OUTPUT_DIR}/vcfs/${PREFIX}_filter6.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${PREFIX}_filter6.vcf.gz"
  echo "INFO: Finished blacklist tagging!"
fi
 
mv "${CURRENT_VCF}" "${OUTPUT_DIR}/vcfs/${PREFIX}.vcf.gz"
mv "${CURRENT_VCF}.tbi" "${OUTPUT_DIR}/vcfs/${PREFIX}.vcf.gz.tbi"

rm -f ${OUTPUT_DIR}/vcfs/${PREFIX}_*chr* \
      ${OUTPUT_DIR}/vcfs/${PREFIX}_raw.vcf.gz* \
      ${OUTPUT_DIR}/vcfs/${PREFIX}_reordered.vcf.gz* \
      ${OUTPUT_DIR}/vcfs/${PREFIX}_filter?.vcf.gz* \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_*chr* \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_blacklist_header.txt

echo "SUCCESS"