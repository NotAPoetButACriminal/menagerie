#!/bin/bash
#
#SBATCH -J cohorc
#SBATCH --nodes 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 500G
#SBATCH --time 3-00:00:00

# --- Usage function ---
usage() {
  cat <<EOF
This script creates a filtered, joint-genotyped cohort VCF file from per-sample GVCF files.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] cohorc.sh -I <gvcfs> -O <out_dir> -C <cohort> -R <ref.fa> [-L <intervals.bed>] ...

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
By default the script will use 64 threads to utilize per chromosome parallelization and speed up joint genotyping.
Do not change the number of threads unless using intervals (-L).
The best way to name log files is cohort name followed by job name and ID (-o  .../logs/COHORT_%x_%A.log).

Variant filtering is done with VQSR, which needs a reasonably large cohort to train on.
Cohorts of fewer than 10 samples are rejected. Trio calling and pedigree-based genotype
refinement are out of scope for this script.

Mandatory flags:
  -I <input>     The GVCF files to joint-genotype. You can use the output of varwolf.sh --gvcf.
                 Either a sample sheet with one .g.vcf.gz path per line.
                 Example sample sheet content:
                   /path/to/SAMPLE1.g.vcf.gz
                   /path/to/SAMPLE2.g.vcf.gz
                 Or the same paths given directly on the command line, separated by commas.
                 Example:
                   /path/to/SAMPLE1.g.vcf.gz,/path/to/SAMPLE2.g.vcf.gz
  -O <dir>       Path to the desired output directory.
  -C <name>      The name of the cohort, used for naming output files.
  -R <file>      Path to the reference genome FASTA file the GVCFs were called against.
                 Must have a .fai index and a .dict sequence dictionary alongside it.

Optional flags:
  --build <name>         Resource set to use for VQSR. Currently only 'hg38' (default).
                         The reference FASTA is NOT set by this flag, since different hg38
                         subversions (with or without alt contigs) all work against the same resources.
  -L <file>              Path to a BED file with target intervals for exome or targeted sequencing.
                         This disables per chromosome parallelization so make sure to give the job
                         2 threads (sbatch -c 2) to not waste cpu.
  --snp-ts <float>       SNP truth sensitivity filter level for ApplyVQSR (default 99.5).
  --indel-ts <float>     INDEL truth sensitivity filter level for ApplyVQSR (default 99.0).
  --max-gaussians <int>  Maximum Gaussians for the VQSR INDEL model (default 4).
                         Lower this if VariantRecalibrator fails to converge on a small cohort.
  --excess-het <float>   ExcessHet phred threshold above which sites are filtered (default 54.69).
  --no-excess-het        Skip the ExcessHet filter entirely. The default threshold assumes
                         unrelated samples and is not appropriate for cohorts containing relatives.
  --min-dp <int>         Genotypes with depth below this are set to no-call (default 10).
  --min-gq <int>         Genotypes with genotype quality below this are set to no-call (default 20).
  --min-ab <float>       Heterozygous genotypes with allelic balance below this are set to no-call (default 0.2).
  --no-genotype-filter   Skip sample level genotype filtering entirely.
  --max-missing <float>  Drop variants missing in more than this fraction of samples (default 0.25).
  --no-normalize         Skip the normalization stage (left align, split multiallelics,
                         deduplicate, keep biallelics, drop missing and monomorphic variants).
  --update-gdb           Add the input samples to existing GenomicsDB workspaces instead of
                         erroring out. Without this, an existing workspace is never touched.
  --keep-intermediates   Keep per chromosome shards and intermediate VCFs instead of deleting them.

Output:
  <cohort>.vcf.gz            All sites, with ExcessHet and VQSR results recorded in the FILTER column.
  <cohort>_<filters>.vcf.gz  PASS sites only, genotype filtered and normalized. The suffix lists
                             the thresholds that were actually applied.
  GenomicsDB workspaces under <out_dir>/gdbs/<cohort>/ are always kept so that --update-gdb
  can add samples later without rebuilding them.

EOF
  exit 1
}

# --- Initial check ---
if [ "$#" -eq 0 ]; then
  usage
fi

set -euo pipefail

# --- Argument Parsing ---
INPUT_GVCFS=""
OUTPUT_DIR=""
COHORT=""
REF=""
BUILD="hg38"
INTERVAL_FILE=""
INTERVALS=""
SNP_TS="99.5"
INDEL_TS="99.0"
MAX_GAUSSIANS="4"
EXCESS_HET="54.69"
RUN_EXCESS_HET=true
MIN_DP="10"
MIN_GQ="20"
MIN_AB="0.2"
RUN_GENOTYPE_FILTER=true
MAX_MISSING="0.25"
RUN_NORMALIZE=true
UPDATE_GDB=false
KEEP_INTERMEDIATES=false

# Manual loop to process options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -I) INPUT_GVCFS="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -C) COHORT="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    --build) BUILD="$2"; shift 2 ;;
    -L)
      INTERVAL_FILE="$2"
      if [[ "${INTERVAL_FILE}" != *.bed ]]; then echo "Error: -L file must be a .bed file." >&2; usage; fi
      if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: ${INTERVAL_FILE}" >&2; usage; fi
      INTERVALS="-L ${INTERVAL_FILE}"
      shift 2
      ;;
    --snp-ts) SNP_TS="$2"; shift 2 ;;
    --indel-ts) INDEL_TS="$2"; shift 2 ;;
    --max-gaussians) MAX_GAUSSIANS="$2"; shift 2 ;;
    --excess-het) EXCESS_HET="$2"; shift 2 ;;
    --no-excess-het) RUN_EXCESS_HET=false; shift ;;
    --min-dp) MIN_DP="$2"; shift 2 ;;
    --min-gq) MIN_GQ="$2"; shift 2 ;;
    --min-ab) MIN_AB="$2"; shift 2 ;;
    --no-genotype-filter) RUN_GENOTYPE_FILTER=false; shift ;;
    --max-missing) MAX_MISSING="$2"; shift 2 ;;
    --no-normalize) RUN_NORMALIZE=false; shift ;;
    --update-gdb) UPDATE_GDB=true; shift ;;
    --keep-intermediates) KEEP_INTERMEDIATES=true; shift ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$INPUT_GVCFS" ]]; then echo "Error: -I <gvcfs> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$COHORT" ]]; then echo "Error: -C <cohort_name> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi

# Validate the reference and its companion index files
if [[ ! -f "$REF" ]]; then echo "Error: Reference FASTA not found: ${REF}" >&2; exit 1; fi
if [[ ! -f "${REF}.fai" ]]; then
  echo "Error: Reference index not found: ${REF}.fai. Create it with 'samtools faidx ${REF}'." >&2
  exit 1
fi
if [[ ! -f "${REF%.*}.dict" ]]; then
  echo "Error: Sequence dictionary not found: ${REF%.*}.dict. Create it with 'gatk CreateSequenceDictionary -R ${REF}'." >&2
  exit 1
fi

# Validate numeric arguments
for ARG_PAIR in "--snp-ts:${SNP_TS}" "--indel-ts:${INDEL_TS}" "--excess-het:${EXCESS_HET}" "--min-ab:${MIN_AB}" "--max-missing:${MAX_MISSING}"; do
  ARG_NAME="${ARG_PAIR%%:*}"
  ARG_VALUE="${ARG_PAIR#*:}"
  if ! [[ "$ARG_VALUE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: ${ARG_NAME} must be a non-negative number (e.g. 99.5)." >&2
    exit 1
  fi
done
for ARG_PAIR in "--max-gaussians:${MAX_GAUSSIANS}" "--min-dp:${MIN_DP}" "--min-gq:${MIN_GQ}"; do
  ARG_NAME="${ARG_PAIR%%:*}"
  ARG_VALUE="${ARG_PAIR#*:}"
  if ! [[ "$ARG_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: ${ARG_NAME} must be a non-negative integer (e.g. 10)." >&2
    exit 1
  fi
done
if ! [[ "$MIN_AB" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]]; then
  echo "Error: --min-ab must be a number between 0 and 1 (e.g. 0.2 for 20%)." >&2
  exit 1
fi
if ! [[ "$MAX_MISSING" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]]; then
  echo "Error: --max-missing must be a number between 0 and 1 (e.g. 0.25 for 25%)." >&2
  exit 1
fi

# --- Resource sets ---
# Each build carries both its VQSR resource files and its chromosome naming convention.
case "$BUILD" in
  hg38)
    HAPMAP="/lustre/imgge/lab01/refs/db/hg38/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
    OMNI="/lustre/imgge/lab01/refs/db/hg38/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
    ONEKG="/lustre/imgge/lab01/refs/db/hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    MILLS="/lustre/imgge/lab01/refs/db/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"
    CHRS=(chr{1..22} chrX chrY chrM)
    ;;
  *)
    echo "Error: --build must be 'hg38'. Got '${BUILD}'." >&2
    usage
    ;;
esac

for RESOURCE in "$HAPMAP" "$OMNI" "$ONEKG" "$MILLS" "$DBSNP"; do
  if [[ ! -f "$RESOURCE" ]]; then
    echo "Error: ${BUILD} VQSR resource not found: ${RESOURCE}" >&2
    exit 1
  fi
done

# --- Start script ---
# Create output directories
mkdir -p "${OUTPUT_DIR}/vcfs/metrics/${COHORT}"
mkdir -p "${OUTPUT_DIR}/gdbs"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

# --- Resolve input GVCFs ---
# -I accepts either a sample sheet file or a comma separated list on the command line.
if [ -f "$INPUT_GVCFS" ]; then
  echo "INFO: Reading GVCF paths from sample sheet ${INPUT_GVCFS}."
  mapfile -t GVCF_PATHS < <(grep -v '^[[:space:]]*$' "$INPUT_GVCFS" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
else
  echo "INFO: Reading GVCF paths from the command line."
  IFS=',' read -r -a GVCF_PATHS <<< "$INPUT_GVCFS"
fi

for GVCF in "${GVCF_PATHS[@]}"; do
  if [[ ! -f "$GVCF" ]]; then echo "Error: GVCF file not found: ${GVCF}" >&2; exit 1; fi
  if [[ ! -f "${GVCF}.tbi" ]]; then
    echo "Error: GVCF index not found: ${GVCF}.tbi. Create it with 'tabix ${GVCF}'." >&2
    exit 1
  fi
done

NUM_SAMPLES=${#GVCF_PATHS[@]}
if [ "$NUM_SAMPLES" -lt 10 ]; then
  echo "Error: Only ${NUM_SAMPLES} GVCF(s) provided. VQSR needs a reasonably large cohort to train on," >&2
  echo "       so this script requires at least 10 samples. For smaller cohorts, call each sample" >&2
  echo "       individually with varwolf.sh, which applies hard filters instead." >&2
  exit 1
fi
echo "INFO: Joint genotyping cohort '${COHORT}' from ${NUM_SAMPLES} GVCF files."

# --- Build the sample name map ---
# GenomicsDBImport reads sample names and paths from this file rather than a long -V string,
# which is what GATK recommends for cohorts of this size.
SAMPLE_MAP="${OUTPUT_DIR}/gdbs/${COHORT}_sample_map.tsv"
echo "INFO: Building sample name map..."
: > "${SAMPLE_MAP}"
for GVCF in "${GVCF_PATHS[@]}"; do
  SAMPLE_NAME=$(bcftools query -l "${GVCF}" | head -1)
  if [[ -z "$SAMPLE_NAME" ]]; then
    echo "Error: Could not read a sample name from ${GVCF}." >&2
    exit 1
  fi
  printf '%s\t%s\n' "${SAMPLE_NAME}" "$(readlink -f "${GVCF}")" >> "${SAMPLE_MAP}"
done

DUPLICATES=$(cut -f1 "${SAMPLE_MAP}" | sort | uniq -d)
if [[ -n "$DUPLICATES" ]]; then
  echo "Error: Duplicate sample names found across the input GVCFs:" >&2
  echo "${DUPLICATES}" | sed 's/^/         /' >&2
  exit 1
fi
echo "INFO: Wrote ${NUM_SAMPLES} unique samples to ${SAMPLE_MAP}!"

# --- Joint Genotyping ---
# With -L the whole cohort is imported into a single workspace, otherwise one workspace per
# chromosome is built in parallel.
if [[ -n "$INTERVAL_FILE" ]]; then
  SHARDS=("${COHORT}")
  SHARD_INTERVALS=("${INTERVALS}")
else
  SHARDS=("${CHRS[@]}")
  SHARD_INTERVALS=()
  for CHR in "${CHRS[@]}"; do
    SHARD_INTERVALS+=("-L ${CHR}")
  done
fi

echo "INFO: Importing GVCFs into GenomicsDB..."
for INDEX in "${!SHARDS[@]}"; do
  SHARD="${SHARDS[$INDEX]}"
  SHARD_L="${SHARD_INTERVALS[$INDEX]}"
  GDB_PATH="${OUTPUT_DIR}/gdbs/${COHORT}/${COHORT}_${SHARD}_gdb"

  # An existing workspace is only ever written to when --update-gdb says so.
  if [ -d "${GDB_PATH}" ]; then
    if [ "$UPDATE_GDB" = false ]; then
      echo "Error: GenomicsDB workspace already exists: ${GDB_PATH}" >&2
      echo "       Pass --update-gdb to add these samples to it, or use a different -C or -O." >&2
      exit 1
    fi
    GDB_FLAG="--genomicsdb-update-workspace-path"
  else
    GDB_FLAG="--genomicsdb-workspace-path"
  fi

  (
    gatk GenomicsDBImport \
      ${GDB_FLAG} "${GDB_PATH}" \
      -R "${REF}" \
      --sample-name-map "${SAMPLE_MAP}" \
      ${SHARD_L}
    echo "INFO: Finished GenomicsDBImport for ${SHARD}!"
  ) &
done
wait
echo "INFO: All GenomicsDBImport jobs finished!"

echo "INFO: Genotyping GVCFs..."
SHARD_VCFS=()
for INDEX in "${!SHARDS[@]}"; do
  SHARD="${SHARDS[$INDEX]}"
  SHARD_L="${SHARD_INTERVALS[$INDEX]}"
  SHARD_VCFS+=("-I" "${OUTPUT_DIR}/vcfs/${COHORT}_${SHARD}.vcf.gz")
  (
    gatk GenotypeGVCFs \
      -R "${REF}" \
      -V "gendb://${OUTPUT_DIR}/gdbs/${COHORT}/${COHORT}_${SHARD}_gdb" \
      ${SHARD_L} \
      -O "${OUTPUT_DIR}/vcfs/${COHORT}_${SHARD}.vcf.gz"
    echo "INFO: Finished GenotypeGVCFs for ${SHARD}!"
  ) &
done
wait
echo "INFO: All GenotypeGVCFs jobs finished!"

echo "INFO: Merging shard VCFs..."
gatk MergeVcfs \
  "${SHARD_VCFS[@]}" \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}_raw.vcf.gz"
echo "INFO: Successfully merged VCF!"

CURRENT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}_raw.vcf.gz"

# --- Excess Heterozygosity Filter ---
# Applied before VQSR, matching the GATK joint genotyping workflow.
if [ "$RUN_EXCESS_HET" = true ]; then
  echo "INFO: Tagging sites with ExcessHet > ${EXCESS_HET}..."
  gatk VariantFiltration \
    -V "${CURRENT_VCF}" \
    -filter "ExcessHet > ${EXCESS_HET}" --filter-name "ExcessHet" \
    -O "${OUTPUT_DIR}/vcfs/${COHORT}_excesshet.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}_excesshet.vcf.gz"
  echo "INFO: Finished ExcessHet tagging!"
else
  echo "INFO: Skipping ExcessHet filter (--no-excess-het set)."
fi

# --- Variant Quality Score Recalibration ---
# VQSR only reads site level annotations, so the model is trained on a sites-only VCF.
echo "INFO: Creating sites-only VCF for recalibration..."
gatk MakeSitesOnlyVcf \
  -I "${CURRENT_VCF}" \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}_sitesonly.vcf.gz"
echo "INFO: Created sites-only VCF!"

echo "INFO: Building SNP recalibration model..."
gatk VariantRecalibrator \
  -R "${REF}" \
  -V "${OUTPUT_DIR}/vcfs/${COHORT}_sitesonly.vcf.gz" \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${HAPMAP}" \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 "${OMNI}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${ONEKG}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${DBSNP}" \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}_snp.recal" \
  --tranches-file "${OUTPUT_DIR}/vcfs/${COHORT}_snp.tranches"
echo "INFO: Finished SNP recalibration model!"

# MQ is deliberately left out of the INDEL model, which GATK advises against using it for.
echo "INFO: Building INDEL recalibration model..."
gatk VariantRecalibrator \
  -R "${REF}" \
  -V "${OUTPUT_DIR}/vcfs/${COHORT}_sitesonly.vcf.gz" \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 "${MILLS}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${DBSNP}" \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode INDEL \
  --max-gaussians "${MAX_GAUSSIANS}" \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}_indel.recal" \
  --tranches-file "${OUTPUT_DIR}/vcfs/${COHORT}_indel.tranches"
echo "INFO: Finished INDEL recalibration model!"

echo "INFO: Applying SNP recalibration at ${SNP_TS}% truth sensitivity..."
gatk ApplyVQSR \
  -R "${REF}" \
  -V "${CURRENT_VCF}" \
  --truth-sensitivity-filter-level "${SNP_TS}" \
  --tranches-file "${OUTPUT_DIR}/vcfs/${COHORT}_snp.tranches" \
  --recal-file "${OUTPUT_DIR}/vcfs/${COHORT}_snp.recal" \
  -mode SNP \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}_snp.vcf.gz"
echo "INFO: Finished applying SNP recalibration!"

echo "INFO: Applying INDEL recalibration at ${INDEL_TS}% truth sensitivity..."
gatk ApplyVQSR \
  -R "${REF}" \
  -V "${OUTPUT_DIR}/vcfs/${COHORT}_snp.vcf.gz" \
  --truth-sensitivity-filter-level "${INDEL_TS}" \
  --tranches-file "${OUTPUT_DIR}/vcfs/${COHORT}_indel.tranches" \
  --recal-file "${OUTPUT_DIR}/vcfs/${COHORT}_indel.recal" \
  -mode INDEL \
  -O "${OUTPUT_DIR}/vcfs/${COHORT}.vcf.gz"
echo "INFO: Finished applying INDEL recalibration!"

# First kept output: every site retained, filtering recorded in the FILTER column.
COHORT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}.vcf.gz"

echo "INFO: Collecting metrics for ${COHORT}.vcf.gz..."
bcftools stats \
  -s - \
  "${COHORT_VCF}" \
  > "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}.stats"
# Plotting is optional. It depends on python/matplotlib, which is not always working in the
# conda environment, and a broken plot must not throw away a finished cohort.
plot-vcfstats -P \
  -p "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}_plots" \
  "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}.stats" \
  || echo "WARNING: plot-vcfstats failed. The .stats file is still available." >&2
echo "INFO: Finished collecting metrics!"

# --- Site level filtering ---
# The output name is assembled from the thresholds that were actually applied, so it can never
# drift away from what the file contains.
SUFFIX="VQSRsnp${SNP_TS}_VQSRindel${INDEL_TS}"

echo "INFO: Keeping only PASS sites..."
bcftools view \
  --threads 8 \
  -f PASS \
  "${COHORT_VCF}" \
  -Oz -o "${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
tabix -f "${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
CURRENT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
echo "INFO: Finished keeping PASS sites!"

# --- Sample level filtering ---
# Genotypes that fail depth, quality or allelic balance are converted to no-call.
if [ "$RUN_GENOTYPE_FILTER" = true ]; then
  echo "INFO: Setting genotypes with DP < ${MIN_DP}, GQ < ${MIN_GQ} or allelic balance < ${MIN_AB} to no-call..."
  SUFFIX="${SUFFIX}_gtDP${MIN_DP}_gtGQ${MIN_GQ}_gtAB${MIN_AB}"
  gatk VariantFiltration \
    -V "${CURRENT_VCF}" \
    --genotype-filter-expression "DP < ${MIN_DP}" --genotype-filter-name "LowDP" \
    --genotype-filter-expression "GQ < ${MIN_GQ}" --genotype-filter-name "LowGQ" \
    --genotype-filter-expression "isHet == 1 && AD[0] + AD[1] > 0 && (AD[1].floatValue() / (AD[0].floatValue() + AD[1].floatValue())) < ${MIN_AB}" --genotype-filter-name "LowAB" \
    --set-filtered-genotype-to-no-call true \
    -O "${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
  echo "INFO: Finished genotype filtering!"
else
  echo "INFO: Skipping genotype filtering (--no-genotype-filter set)."
fi

# --- Normalization ---
# Left align and normalize indels, split multiallelics into biallelics, deduplicate, remove "*"
# alleles and keep only biallelics, remove variants missing in too many samples and variants that
# lost their ALT allele.
if [ "$RUN_NORMALIZE" = true ]; then
  echo "INFO: Normalizing and removing variants missing in more than ${MAX_MISSING} of samples..."
  SUFFIX="${SUFFIX}_norm_miss${MAX_MISSING}"
  bcftools norm \
    --threads 8 \
    -f "${REF}" \
    -m-any \
    "${CURRENT_VCF}" \
    -Ou | \
  bcftools norm \
    --threads 8 \
    -f "${REF}" \
    -d all \
    -Ou | \
  bcftools view \
    --threads 8 \
    -e 'ALT="*"' \
    -m2 -M2 \
    -Ou | \
  bcftools filter \
    --threads 8 \
    -e "F_MISSING > ${MAX_MISSING}" \
    -Ou | \
  bcftools +fill-tags \
    -Ou -- -t AC,AN | \
  bcftools view \
    --threads 8 \
    -e 'AC=0 || AC=AN' \
    -Oz -o "${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
  tabix -f "${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
  CURRENT_VCF="${OUTPUT_DIR}/vcfs/${COHORT}_${SUFFIX}.vcf.gz"
  echo "INFO: Finished normalizing!"
else
  echo "INFO: Skipping normalization (--no-normalize set)."
fi

# Second kept output: the analysis ready cohort VCF.
FINAL_VCF="${CURRENT_VCF}"

echo "INFO: Collecting metrics for ${COHORT}_${SUFFIX}.vcf.gz..."
bcftools stats \
  -s - \
  "${FINAL_VCF}" \
  > "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}_${SUFFIX}.stats"
plot-vcfstats -P \
  -p "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}_${SUFFIX}_plots" \
  "${OUTPUT_DIR}/vcfs/metrics/${COHORT}/${COHORT}_${SUFFIX}.stats" \
  || echo "WARNING: plot-vcfstats failed. The .stats file is still available." >&2
echo "INFO: Finished collecting metrics!"

# --- Cleanup ---
# The GenomicsDB workspaces are never removed, so that --update-gdb can reuse them.
if [ "$KEEP_INTERMEDIATES" = false ]; then
  echo "INFO: Removing intermediate files..."
  for FILE in "${OUTPUT_DIR}/vcfs/${COHORT}"_*.vcf.gz "${OUTPUT_DIR}/vcfs/${COHORT}"_*.vcf.gz.tbi; do
    if [[ -f "$FILE" && "$FILE" != "${FINAL_VCF}" && "$FILE" != "${FINAL_VCF}.tbi" ]]; then
      rm -f "$FILE"
    fi
  done
  rm -f "${OUTPUT_DIR}/vcfs/${COHORT}"_*.recal \
        "${OUTPUT_DIR}/vcfs/${COHORT}"_*.recal.idx \
        "${OUTPUT_DIR}/vcfs/${COHORT}"_*.tranches \
        "${OUTPUT_DIR}/vcfs/${COHORT}"_*.tranches.pdf
  echo "INFO: Finished removing intermediate files!"
else
  echo "INFO: Keeping intermediate files (--keep-intermediates set)."
fi

echo "INFO: Cohort VCF (all sites, filters tagged): ${COHORT_VCF}"
echo "INFO: Analysis ready VCF: ${FINAL_VCF}"
echo "INFO: GenomicsDB workspaces kept in: ${OUTPUT_DIR}/gdbs/${COHORT}/"

echo "SUCCESS"
