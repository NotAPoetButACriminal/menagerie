#!/bin/bash
#
#SBATCH -J varwolf
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 128G
#SBATCH --time 3-00:00:00

# --- Usage function ---
usage() {
  cat <<EOF
This script creates a filtered germline VCF file from a BAM file for genome, exome or targeted sequencing.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] varwolf.sh -I <input.bam> -O <out_dir> -S <sample> -R <ref.fa> [-L <intervals.bed>] ...

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
By default the script will use 32 threads to utilize per chromosome parallelization and speed up variant calling.
Do not change the number of threads unless using intervals (-L) or in --singlethread mode.
The best way to name log files is sample name followed by job name and ID (-o  .../logs/SAMPLE_%x_%A.log).

Mandatory flags:
  -I <file>      Path to the input BAM file. You can use the output of bampire.sh.
  -O <dir>       Path to the desired output directory.
  -R <file>      Path to the reference genome FASTA file to which the BAM file is aligned.

Optional flags:
  -S <name>      The name of the sample, used for naming output files. If not provided, it will be derived from the BAM file name.
  -L <file>      Path to a BED file with target intervals for exome or targeted sequencing.
                 This disables multi-threading so make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.
  --gvcf         Also produce a GVCF file that can be used for joint genotyping with cohortcrawler.sh
  --counts       Also produce a read counts HDF5 file for CNV calling with copycat.sh
  --singlethread Disable per chromosome parallelization for genome sequencing and perform variant calling on a single cpu thread.
                 This will be significantly slower but can be used when running a large number of samples in parallel (>50).
                 Make sure to give the job 2 threads (sbatch -c 2) to not waste cpu.

EOF
  exit 1
}

# --- Initial check ---
if [ "$#" -eq 0 ]; then
  usage
fi

set -euo pipefail

# --- Argument Parsing ---
INPUT_BAM=""
OUTPUT_DIR=""
SAMPLE=""
REF=""
INTERVAL_FILE=""
INTERVALS=""
GVCF_MODE=false
SINGLETHREAD_MODE=false
RUN_COLLECT_COUNTS=false

# Manual loop to process options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -I) INPUT_BAM="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    -S) SAMPLE="$2"; shift 2 ;;
    -L)
      INTERVAL_FILE="$2"
      if [[ "${INTERVAL_FILE}" != *.bed ]]; then echo "Error: -L file must be a .bed file." >&2; usage; fi
      if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: ${INTERVAL_FILE}" >&2; usage; fi
      INTERVALS="-L ${INTERVAL_FILE}"
      shift 2
      ;;
    --gvcf) GVCF_MODE=true; shift ;;
    --counts) RUN_COLLECT_COUNTS=true; shift ;;
    --singlethread) SINGLETHREAD_MODE=true; shift ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$INPUT_BAM" ]]; then echo "Error: -I <input.bam> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi
# Derive sample name from BAM file if not provided
if [[ -z "$SAMPLE" ]]; then
  echo "INFO: Sample name not provided with -S. Deriving from input BAM file name."
  SAMPLE=$(basename "${INPUT_BAM}" .bam)
fi

# --- Start script ---
# Create output directories
mkdir -p "${OUTPUT_DIR}/vcfs"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

# Set HaplotypeCaller mode
HC_MODE=""
VCF_SUFFIX=""

if [ "$GVCF_MODE" = true ]; then
  echo "INFO: GVCF mode enabled."
  HC_MODE="-ERC GVCF"
  VCF_SUFFIX=".g.vcf.gz"
else
  VCF_SUFFIX="_raw.vcf.gz"
fi

# --- Germline SNV Calling ---
# Run HaplotypeCaller
if [[ "$SINGLETHREAD_MODE" = true || -n "$INTERVAL_FILE" ]]; then
  echo "INFO: Running HaplotypeCaller in single-thread mode..."
  gatk HaplotypeCaller \
    -R "${REF}" \
    ${INTERVALS} \
    ${HC_MODE} \
    -ip 20 \
    -I "${INPUT_BAM}" \
    -O "${OUTPUT_DIR}/vcfs/${SAMPLE}${VCF_SUFFIX}"
  echo "INFO: Finished HaplotypeCaller in single-thread mode!"
else
  CHRS=(chr{1..22} chrX chrY chrM)
  CHR_VCFS=()
  echo "INFO: Running HaplotypeCaller per chromosome..."
  for CHR in "${CHRS[@]}"; do
    CHR_VCFS+=("-I" "${OUTPUT_DIR}/vcfs/${SAMPLE}_${CHR}${VCF_SUFFIX}")
    (
      gatk HaplotypeCaller \
        -R "${REF}" \
        -L "${CHR}" \
        ${HC_MODE} \
        -I "${INPUT_BAM}" \
        -O "${OUTPUT_DIR}/vcfs/${SAMPLE}_${CHR}${VCF_SUFFIX}"
      echo "INFO: Finished HaplotypeCaller for ${CHR}!"
    ) &
  done
  wait
  echo "INFO: All chromosome jobs finished!"

  echo "INFO: Merging per-chromosome VCFs..."
  gatk MergeVcfs \
    "${CHR_VCFS[@]}" \
    -O "${OUTPUT_DIR}/vcfs/${SAMPLE}${VCF_SUFFIX}"
  echo "INFO:Successfully merged VCF!"
fi

# If in --gvcf mode, create VCF
if [ "$GVCF_MODE" = true ]; then
  echo "INFO: Creating VCF from GVCF..."
  gatk GenotypeGVCFs \
    -R "${REF}" \
    ${INTERVALS} \
    -ip 50 \
    -V "${OUTPUT_DIR}/vcfs/${SAMPLE}${VCF_SUFFIX}" \
    -O "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.vcf.gz"
  echo "INFO: Created VCF from GVCF!"
fi

# --- Variant Filtering ---
echo "INFO: Filtering VCF..."
gatk VariantFiltration \
	  -V "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.vcf.gz" \
    -filter "DP < 5.0" --filter-name "DP5" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${OUTPUT_DIR}/vcfs/${SAMPLE}.vcf.gz"
echo "INFO: Finished filtering VCF!"

# --- Collect Read Counts (Optional) ---

if [ "$RUN_COLLECT_COUNTS" = true ]; then
  mkdir -p "${OUTPUT_DIR}/counts"

  if [[ -n "$INTERVAL_FILE" ]]; then
    echo "INFO: Preprocessing provided intervals for read counting..."
    gatk PreprocessIntervals \
      -R "${REF}" \
      -L "${INTERVAL_FILE}" \
      --bin-length 0 \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O "${OUTPUT_DIR}/counts/${SAMPLE}_bins.interval_list"
    echo "INFO: Created bins!"
  else
    echo "INFO: No interval file provided. Creating whole-genome bins for read counting..."
    gatk PreprocessIntervals \
      -R "${REF}" \
      --bin-length 1000 \
      --padding 0 \
      -imr OVERLAPPING_ONLY \
      -O "${OUTPUT_DIR}/counts/${SAMPLE}_bins.interval_list"
    echo "INFO: Created bins!"
  fi

  echo "INFO: Running CollectReadCounts..."
  gatk CollectReadCounts \
    -R "${REF}" \
    -L "${OUTPUT_DIR}/counts/${SAMPLE}_bins.interval_list" \
    -I "${INPUT_BAM}" \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O "${OUTPUT_DIR}/counts/${SAMPLE}.hdf5"
  echo "INFO: Finished collecting read counts!"
fi

rm ${OUTPUT_DIR}/vcfs/${SAMPLE}_*

echo "SUCCESS"
