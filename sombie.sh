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
This script creates a filtered somatic VCF file from matched tumor/normal BAM files using GATK Mutect2.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] sombie.sh -T <tumor.bam> -N <normal.bam> -O <out_dir> -R <ref.fa> [-L <intervals.bed>] ...

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
By default the script will use 32 threads to utilize per chromosome parallelization and speed up variant calling.
Do not change the number of threads unless using intervals (-L) or in --singlethread mode.
The best way to name log files is sample name followed by job name and ID (-o  .../logs/SAMPLE_%x_%A.log).

Mandatory flags:
  -T <file>      Path to the tumor BAM file. You can use the output of bampire.sh.
  -N <file>      Path to the matched normal BAM file. You can use the output of bampire.sh.
  -O <dir>       Path to the desired output directory.
  -R <file>      Path to the reference genome FASTA file to which the BAM files are aligned.

Optional flags:
  -S <name>      Tumor sample name, used for naming output files. If not provided, it is read from the tumor BAM's read group (SM tag).
  -NS <name>     Normal sample name, passed to Mutect2 as the -normal argument. If not provided, it is read from the normal BAM's read group (SM tag).
  -L <file>      Path to a BED file with target intervals for exome or targeted sequencing.
                 This disables multi-threading so make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.
  --custom-pon <file>        Path to a custom Panel of Normals VCF. Defaults to the GATK-provided 1000 Genomes PoN
                              (see DEFAULT_PON near the top of this script). Point this at your own cohort-built PoN,
                              or a merged/unioned PoN, if you have one.
  --custom-germline <file>   Path to a custom germline resource VCF (e.g. gnomAD af-only). Defaults to the GATK-provided
                              gnomAD resource (see DEFAULT_GERMLINE near the top of this script).
  --custom-common <file>     Path to a custom small biallelic common-SNP VCF for contamination estimation. Defaults to
                              the GATK-provided gnomAD biallelic subset (see DEFAULT_COMMON near the top of this script).
  --genotype-pon-sites       Emit and genotype sites that match the panel of normals instead of dropping them, tagging
                              them with a PoN filter so they remain visible (but excluded from PASS) in the output VCF.
  --genotype-germline-sites  Emit and genotype sites that look germline per the germline resource instead of dropping
                              them, tagging them with a germline filter so they remain visible (but excluded from PASS).
  --normal-lod <float>       Override Mutect2's --normal-lod threshold (default 2.2). Higher values make it harder to
                              reject a candidate based on normal-sample evidence (more permissive toward calling somatic).
  --disable-mate-filter      Disable MateOnSameContigOrNoMappedMateReadFilter, retaining read pairs whose mate maps to
                              a different contig. Useful for capturing evidence near translocation breakpoints; adds noise.
  --singlethread Disable per chromosome parallelization for genome sequencing and perform variant calling on a single cpu thread.
                 This will be significantly slower but can be used when running a large number of samples in parallel (>50).
                 Make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.

EOF
  exit 1
}

# --- Default resources ---
DEFAULT_PON="/lustre/imgge/lab01/refs/db/hg38/gatk_1000g_pon.hg38.vcf.gz"
DEFAULT_GERMLINE="/lustre/imgge/lab01/refs/db/hg38/gatk_af-only-gnomad.hg38.vcf.gz"
DEFAULT_COMMON="/lustre/imgge/lab01/refs/db/hg38/gatk_af-only-gnomad-common-biallelic.vcf.gz"

# --- Initial check ---
if [ "$#" -eq 0 ]; then
  usage
fi

set -euo pipefail

# --- Argument Parsing ---
TUMOR_BAM=""
NORMAL_BAM=""
OUTPUT_DIR=""
REF=""
PON="${DEFAULT_PON}"
GERMLINE_RESOURCE="${DEFAULT_GERMLINE}"
CONTAM_RESOURCE="${DEFAULT_COMMON}"
TUMOR_SAMPLE=""
NORMAL_SAMPLE=""
INTERVAL_FILE=""
INTERVALS=""
GENOTYPE_PON=false
GENOTYPE_GERMLINE=false
NORMAL_LOD=""
DISABLE_MATE_FILTER=false
SINGLETHREAD_MODE=false

# Manual loop to process options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -T) TUMOR_BAM="$2"; shift 2 ;;
    -N) NORMAL_BAM="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    --custom-germline) GERMLINE_RESOURCE="$2"; shift 2 ;;
    --custom-pon) PON="$2"; shift 2 ;;
    --custom-common) CONTAM_RESOURCE="$2"; shift 2 ;;
    -S) TUMOR_SAMPLE="$2"; shift 2 ;;
    -NS) NORMAL_SAMPLE="$2"; shift 2 ;;
    -L)
      INTERVAL_FILE="$2"
      if [[ "${INTERVAL_FILE}" != *.bed ]]; then echo "Error: -L file must be a .bed file." >&2; usage; fi
      if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: ${INTERVAL_FILE}" >&2; usage; fi
      INTERVALS="-L ${INTERVAL_FILE} -ip 15"
      shift 2
      ;;
    --genotype-pon-sites) GENOTYPE_PON=true; shift ;;
    --genotype-germline-sites) GENOTYPE_GERMLINE=true; shift ;;
    --normal-lod) NORMAL_LOD="$2"; shift 2 ;;
    --disable-mate-filter) DISABLE_MATE_FILTER=true; shift ;;
    --singlethread) SINGLETHREAD_MODE=true; shift ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$TUMOR_BAM" ]]; then echo "Error: -T <tumor.bam> is a mandatory flag." >&2; usage; fi
if [[ -z "$NORMAL_BAM" ]]; then echo "Error: -N <normal.bam> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi
if [[ ! -f "$GERMLINE_RESOURCE" ]]; then
  echo "Error: Germline resource not found at ${GERMLINE_RESOURCE}" >&2
  echo "       Either provide --custom-germline <file>, or update DEFAULT_GERMLINE near the top of this script to a valid path." >&2
  exit 1
fi
if [[ ! -f "$PON" ]]; then
  echo "Error: Panel of Normals not found at ${PON}" >&2
  echo "       Either provide --custom-pon <file>, or update DEFAULT_PON near the top of this script to a valid path." >&2
  exit 1
fi
if [[ ! -f "$CONTAM_RESOURCE" ]]; then
  echo "Error: Common-SNP contamination resource not found at ${CONTAM_RESOURCE}" >&2
  echo "       Either provide --custom-common <file>, or update DEFAULT_COMMON near the top of this script to a valid path." >&2
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
if [[ -z "$NORMAL_SAMPLE" ]]; then
  echo "INFO: Normal sample name not provided with -NS. Deriving from normal BAM read group (SM tag)."
  NORMAL_SAMPLE=$(samtools view -H "${NORMAL_BAM}" | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/' | head -1)
  if [[ -z "$NORMAL_SAMPLE" ]]; then echo "Error: Could not derive normal sample name from BAM read group." >&2; exit 1; fi
fi
echo "INFO: Tumor sample: ${TUMOR_SAMPLE} | Normal sample: ${NORMAL_SAMPLE}"

PREFIX="${TUMOR_SAMPLE}"

# Build optional Mutect2 argument string
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
    -I "${NORMAL_BAM}" \
    -tumor "${TUMOR_SAMPLE}" \
    -normal "${NORMAL_SAMPLE}" \
    --germline-resource "${GERMLINE_RESOURCE}" \
    --panel-of-normals "${PON}" \
    --f1r2-tar-gz "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz" \
    "${MUTECT_EXTRA_ARGS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz"
  echo "INFO: Finished Mutect2 in single-thread mode!"

  STATS_ARGS=("-stats" "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz.stats")
  F1R2_ARGS=("-I" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz")
  cp "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz" "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw_merged.vcf.gz"
  cp "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz.tbi" "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw_merged.vcf.gz.tbi"
else
  CHRS=(chr{1..22} chrX chrY chrM)
  CHR_VCFS=()
  STATS_ARGS=()
  F1R2_ARGS=()
  echo "INFO: Running Mutect2 per chromosome..."
  for CHR in "${CHRS[@]}"; do
    CHR_VCFS+=("-I" "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}_somatic_raw.vcf.gz")
    STATS_ARGS+=("-stats" "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}_somatic_raw.vcf.gz.stats")
    F1R2_ARGS+=("-I" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_${CHR}_f1r2.tar.gz")
    (
      gatk Mutect2 \
        -R "${REF}" \
        -L "${CHR}" \
        -I "${TUMOR_BAM}" \
        -I "${NORMAL_BAM}" \
        -normal "${NORMAL_SAMPLE}" \
        --germline-resource "${GERMLINE_RESOURCE}" \
        --panel-of-normals "${PON}" \
        --f1r2-tar-gz "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_${CHR}_f1r2.tar.gz" \
        "${MUTECT_EXTRA_ARGS[@]}" \
        -O "${OUTPUT_DIR}/vcfs/${PREFIX}_${CHR}_somatic_raw.vcf.gz"
      echo "INFO: Finished Mutect2 for ${CHR}!"
    ) &
  done
  wait
  echo "INFO: All chromosome jobs finished!"

  echo "INFO: Merging per-chromosome VCFs..."
  gatk MergeVcfs \
    "${CHR_VCFS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw_merged.vcf.gz"
  echo "INFO: Successfully merged VCF!"

  echo "INFO: Merging per-chromosome stats files..."
  gatk MergeMutectStats \
    "${STATS_ARGS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz.stats"
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

echo "INFO: Running GetPileupSummaries for normal..."
gatk GetPileupSummaries \
  -I "${NORMAL_BAM}" \
  -V "${CONTAM_RESOURCE}" \
  -L "${CONTAM_RESOURCE}" \
  ${INTERVALS} \
  -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table"

echo "INFO: Calculating contamination..."
gatk CalculateContamination \
  -I "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table" \
  -matched "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table" \
  -O "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_contamination.table" \
  --tumor-segmentation "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_segments.table"
echo "INFO: Finished contamination estimation!"

FILTER_EXTRA_ARGS+=("--contamination-table" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_contamination.table")
FILTER_EXTRA_ARGS+=("--tumor-segmentation" "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_segments.table")

# --- Variant Filtering ---
echo "INFO: Filtering VCF..."
gatk FilterMutectCalls \
  -R "${REF}" \
  -V "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw_merged.vcf.gz" \
  --ob-priors "${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_read-orientation-model.tar.gz" \
  --stats "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz.stats" \
  "${FILTER_EXTRA_ARGS[@]}" \
  -O "${OUTPUT_DIR}/vcfs/${PREFIX}_somatic.vcf.gz"
echo "INFO: Finished filtering VCF!"

rm -f ${OUTPUT_DIR}/vcfs/${PREFIX}_*chr* \
      ${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw_merged.vcf.gz* \
      ${OUTPUT_DIR}/vcfs/${PREFIX}_somatic_raw.vcf.gz* \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_f1r2.tar.gz \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_tumor_pileups.table \
      ${OUTPUT_DIR}/vcfs/metrics/${PREFIX}_normal_pileups.table

echo "SUCCESS"