#!/bin/bash
#
#SBATCH -J bampire
#SBATCH --nodes 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 256G
#SBATCH --time 1-00:00:00

# --- Usage function ---
usage() {
  cat <<EOF
This script creates an analysis-ready BAM file from raw reads (FASTQ) for genome, exome or targeted sequencing.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] bampire.sh -I <fastq.csv> -O <out_dir> -S <sample> -R <ref.fa>

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
A minimum of 8 threads is required (-c 8). The best way to name log files is sample name followed by job name and ID (-o  .../logs/SAMPLE_%x_%A.log).

Mandatory arguments:
  -I <file>      Path to the pair of FASTQ files, separated by a comma.
                 Example:
                   /path/to/SAMPLE_R1.fastq.gz,/path/to/SAMPLE_R2.fastq.gz
                 If the sample has multiple pairs of FASTQ files split across lanes, provide a CSV file with one pair per line.
                 Example CSV content:
                   /path/to/SAMPLE_R1_L001.fastq.gz,/path/to/SAMPLE_R2_L001.fastq.gz
                   /path/to/SAMPLE_R1_L002.fastq.gz,/path/to/SAMPLE_R2_L002.fastq.gz
  -O <dir>       Path to the desired output directory.
  -S <name>      The name of the sample, used for naming output files.
  -R <file>      Path to the reference genome FASTA file (bwa-mem2 index must be present).

Optional arguments:
  --skip-bqsr    Skip the Base Quality Score Recalibration (BQSR) step. This speeds up run time significantly and decreases file size with minimal difference in output, but it is still not considered best practice by GATK.
  --legacy-bwa   Use the original bwa aligner instead of the default bwa-mem2. Useful if no bwa-mem2 index is available.
  --custom-rg <string>   If the reads were not sequenced with an Illumina sequencer, you will need to provide custom read group information.
                 Example: --custom-rg "@RG\tID:FLOWCELL\tPL:ILLUMINA\tLB:LIBRARY\tSM:SAMPLE"

EOF
  exit 1
}

# --- Initial check ---
# If the script is called with no arguments, display the usage message and exit.
if [ "$#" -eq 0 ]; then
  usage
fi

set -eux

# --- Argument Parsing ---
# Flags
INPUT_FILES=""
OUTPUT_DIR=""
SAMPLE=""
REF=""
CUSTOM_RG_STRING=""
SKIP_BQSR=false
LEGACY_BWA=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -I) INPUT_FILES="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -S) SAMPLE="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    --skip-bqsr) SKIP_BQSR=true; shift ;;
    --legacy-bwa) LEGACY_BWA=true; shift ;;
    --custom-rg) CUSTOM_RG_STRING="$2"; shift 2 ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$INPUT_FILES" ]]; then echo "Error: -I <fastq_list.csv> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$SAMPLE" ]]; then echo "Error: -S <sample_name> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi


# Set up remaining variables
THREADS=${SLURM_CPUS_PER_TASK:-8}
DBSNP="/lustre/imgge/db/hg38/hg38.dbsnp155.vcf.gz"

# --- Start script ---
# Create output directories
mkdir -p "${OUTPUT_DIR}/bams/metrics"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

# Set aligner
if [ "$LEGACY_BWA" = true ]; then
  BWA="bwa"
  echo "Using legacy bwa aligner."
else
  BWA="bwa-mem2"
fi

# Load input files
FASTQ_LIST=""
if [ -f "$INPUT_FILES" ]; then
    FASTQ_LIST=$(cat "$INPUT_FILES")
else
    FASTQ_LIST="$INPUT_FILES"
fi

# Issue warning if sample name is not present in fastq file names
while IFS=',' read -r R1_FILE R2_FILE; do
  if [[ ! "$R1_FILE" == *"$SAMPLE"* ]]; then
    echo "WARNING: Sample name '${SAMPLE}' not found in FASTQ filename '${R1_FILE}'. Please verify your inputs."
  fi
done <<< "$FASTQ_LIST"

# --- Alignment and Preprocessing ---
# Read the CSV file line by line. Each line is one lane (one pair of FASTQ files).
LANE_NUM=0
BAMS_PER_LANE=()
while IFS=',' read -r R1_FILE R2_FILE; do
  LANE_NUM=$((LANE_NUM + 1))
  if [[ -n "$R1_FILE" && -n "$R2_FILE" ]]; then
    RG_STRING=""
    if [[ -n "$CUSTOM_RG_STRING" ]]; then
      RG_STRING="${CUSTOM_RG_STRING}"
    else
      FLOWCELL=$(zcat -f "${R1_FILE}" | head -n 1 | cut -d ":" -f 3)
      LIBRARY=$(zcat -f "${R1_FILE}" | head -n 1 | cut -d ":" -f 2 | sed 's/^/Lib/g')
      RG_STRING="@RG\tID:${FLOWCELL}.${LANE_NUM}\tPL:ILLUMINA\tLB:${LIBRARY}\tSM:${SAMPLE}"
    fi
    echo "Aligning ${SAMPLE} Lane ${LANE_NUM}"
    fastp \
      -w $((THREADS / 8)) \
      -i "${R1_FILE}" \
      -I "${R2_FILE}" \
      --stdout \
      -j "${OUTPUT_DIR}/bams/metrics/${SAMPLE}_fastp.json" \
      -h "${OUTPUT_DIR}/bams/metrics/${SAMPLE}_fastp.html" \
    | ${BWA} mem \
      -t "${THREADS}" \
      -M -p \
      -R "${RG_STRING}" \
      "${REF}" - \
    | samtools sort \
      -n \
      -@ $((THREADS / 8)) \
    > "${OUTPUT_DIR}/bams/${SAMPLE}_L${LANE_NUM}.bam"
    BAMS_PER_LANE+=("-I" "${OUTPUT_DIR}/bams/${SAMPLE}_L${LANE_NUM}.bam")
  fi
done <<< "$FASTQ_LIST"

echo "Finished mapping all reads for ${SAMPLE}!"

echo "Marking duplicates for ${SAMPLE}..."
gatk MarkDuplicatesSpark \
  --java-options "-Xmx32G" \
  -R "${REF}" \
  "${BAMS_PER_LANE[@]}" \
  -O "${OUTPUT_DIR}/bams/${SAMPLE}_dd.bam" \
  -M "${OUTPUT_DIR}/bams/metrics/${SAMPLE}_mdmetrics.txt" \
  --spark-runner LOCAL --spark-master local[${THREADS}]
echo "Finished marking duplicates for ${SAMPLE}!"

# Conditionally run or skip BQSR based on the --skip-bqsr flag.
if [ "$SKIP_BQSR" = true ]; then
  echo "Skipping BQSR step as requested."
  mv "${OUTPUT_DIR}/bams/${SAMPLE}_dd.bam" "${OUTPUT_DIR}/bams/${SAMPLE}.bam"
  mv "${OUTPUT_DIR}/bams/${SAMPLE}_dd.bam.bai" "${OUTPUT_DIR}/bams/${SAMPLE}.bam.bai"
else
  echo "Recalibrating bases for ${SAMPLE}..."
  gatk BQSRPipelineSpark \
    --java-options "-Xmx32G" \
    -R "${REF}" \
    -I "${OUTPUT_DIR}/bams/${SAMPLE}_dd.bam" \
    -O "${OUTPUT_DIR}/bams/${SAMPLE}.bam" \
    --known-sites ${DBSNP} \
    --spark-runner LOCAL --spark-master local[${THREADS}]
  echo "Finished recalibrating bases for ${SAMPLE}!"
fi

rm \
  "${OUTPUT_DIR}/bams/${SAMPLE}"_dd.bam* \
  "${OUTPUT_DIR}/bams/${SAMPLE}"_L*.bam*

echo "SUCCESS"