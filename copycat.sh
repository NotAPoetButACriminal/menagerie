#!/bin/bash
#
#SBATCH -J copycat
#SBATCH --nodes 1
#SBATCH --cpus-per-task 128
#SBATCH --mem 256G
#SBATCH --time 3-00:00:00

# --- Usage function ---
usage() {
  cat <<EOF
This script calls germline copy number variants from read count HDF5 files using GATK gCNV.

Usage: sbatch [-c <num_cpus>] [-o <logfile_path>] copycat.sh -I <hdf5s> -O <out_dir> -R <ref.fa> {-C <cohort> | -M <model_dir>} [-L <intervals.bed>] ...

This script is designed to be run with SLURM using sbatch. DO NOT run it on the login node except for testing.
By default the script will use 128 threads, the interval shards run concurrently and each one multithreads internally.
The best way to name log files is cohort name followed by job name and ID (-o  .../logs/COHORT_%x_%A.log).

The script runs in one of two modes:

  COHORT mode (default) fits a new CNV model to the cohort and calls CNVs in the same samples.
  This is what you want the first time you process a batch of samples.
  The samples must be sequenced and processed the same way.
  GATK needs a reasonably large cohort to fit a usable model, so fewer than 10 samples is rejected.

  CASE mode ( when -M <dir> is provided) calls CNVs in new samples against a model built by an earlier cohort run.
  This allows CNVs to be called for single samples, as long as they were sequenced and processed the same way as the cohort.
  Provide the <out_dir>/gcnv/<cohort> directory that the earlier run produced.

All samples must be counted against an identical interval list. The simplest way to guarantee that
is to run varwolf.sh --counts with the same -L for every sample, then pass that same -L here.

Mandatory flags:
  -I <input>     The read count HDF5 files to call CNVs from. You can use the output of
                 varwolf.sh --counts.
                 In COHORT mode it is either a sample sheet with one .hdf5 path per line.
                 Example sample sheet content:
                   /path/to/SAMPLE1.hdf5
                   /path/to/SAMPLE2.hdf5
                 Or the same paths given directly on the command line, separated by commas.
                 Example:
                   /path/to/SAMPLE1.hdf5,/path/to/SAMPLE2.hdf5
                 In CASE mode this is a single .hdf5 file, since one job calls one sample.
  -O <dir>       Path to the desired output directory.
  -R <file>      Path to the reference genome FASTA file the BAMs were aligned to.
                 Must have a .fai index and a .dict sequence dictionary alongside it.

Cohort mode:
  -C <name>      The name of the cohort, used for naming the model and working directories.
                 Giving it selects COHORT mode, which is the default.

Case mode:
  -M <dir>       The <out_dir>/gcnv/<cohort> directory of an earlier COHORT run for calling in CASE mode.
                 Only samples sequenced and processed the exact same way as a previous cohort
                 can be run in case mode against that cohort's model.
                 Giving it selects CASE mode, and the work directory is named after the sample, taken from
                 the .hdf5 file name. In this mode -L, --scatters and --low-count-pct are ignored, since the
                 intervals and the shard layout are taken from the model.

Optional flags:
  -L <file>              The same BED file that was passed to varwolf.sh -L when the counts were collected.
                         The bins are rebuilt from it exactly the way varwolf.sh built them, so any other
                         file gives intervals that do not line up with the counts.
                         If omitted, whole genome 1000 bp bins are used, again matching varwolf.sh.
  --build <name>         Resource set to use. Currently only 'hg38' (default). This selects the
                         pseudoautosomal regions and the contig ploidy priors, not the reference
                         FASTA, since the various hg38 subversions share both.
  --scatters <int>       Number of shards to split the intervals into for GermlineCNVCaller parallelism
                         (default 10). Use --scatters 1 to skip scattering and call every interval in one job.
  --ploidy-priors <file> Custom contig ploidy priors table for DetermineGermlineContigPloidy.
  --custom-par <file>    BED file of pseudoautosomal regions to exclude, overriding the build default.
  --low-count-pct <int>  FilterIntervals drops intervals with a low count in more than this
                         percentage of samples (default 65). COHORT mode only.
  --min-qual <float>     Segments below this QUAL are tagged CNVQUAL in the FILTER column (default 100).
  --rmv-qual <float>     Segments below this QUAL are dropped from the final VCF (default 30).
  --keep-intermediates   Keep the raw and filter tagged CNV VCFs instead of deleting them.
  --overwrite            Delete the ploidy model, interval shards and gCNV shards of an earlier run with the
                         same name before starting. Without this, a run whose working directory already holds
                         results is refused, because shards from two runs mix and silently corrupt the calls.
                         Overwriting a COHORT model also invalidates any CASE runs called against it.

Output:
  <out_dir>/vcfs/<sample>.cnv.vcf.gz            Final per sample CNV calls, indexed, with SVTYPE=CNV
                                                filled in so downstream SV tools accept them.
  <out_dir>/gcnv/<cohort>/                      Ploidy and gCNV models, scattered calls, denoised copy
                                                ratios and per interval genotypes (<sample>_intervals.cnv.vcf.gz).
                                                Never deleted unless --overwrite is given, so it can be
                                                reused as the -M directory for later CASE runs.
                                                A CASE run writes the same files into <out_dir>/gcnv/<sample>/.

EOF
  exit 1
}

# --- Initial check ---
if [ "$#" -eq 0 ]; then
  usage
fi

set -euo pipefail

# --- Argument Parsing ---
INPUT_COUNTS=""
OUTPUT_DIR=""
COHORT=""
REF=""
MODEL_DIR=""
GENOME_BUILD="hg38"
INTERVAL_FILE=""
PLOIDY_PRIORS=""
CUSTOM_PAR=""
LOW_COUNT_PCT="65"
SCATTERS_REQUESTED="10"
MIN_QUAL="100.0"
RMV_QUAL="30.0"
KEEP_INTERMEDIATES=false
OVERWRITE=false

# Manual loop to process options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -I) INPUT_COUNTS="$2"; shift 2 ;;
    -O) OUTPUT_DIR="$2"; shift 2 ;;
    -C) COHORT="$2"; shift 2 ;;
    -R) REF="$2"; shift 2 ;;
    -M) MODEL_DIR="$2"; shift 2 ;;
    --build) GENOME_BUILD="$2"; shift 2 ;;
    -L)
      INTERVAL_FILE="$2"
      if [[ "${INTERVAL_FILE}" != *.bed ]]; then
        echo "Error: -L file must be a .bed file, the same one that was given to varwolf.sh -L." >&2
        usage
      fi
      if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: ${INTERVAL_FILE}" >&2; usage; fi
      shift 2
      ;;
    --scatters) SCATTERS_REQUESTED="$2"; shift 2 ;;
    --ploidy-priors) PLOIDY_PRIORS="$2"; shift 2 ;;
    --custom-par) CUSTOM_PAR="$2"; shift 2 ;;
    --low-count-pct) LOW_COUNT_PCT="$2"; shift 2 ;;
    --min-qual) MIN_QUAL="$2"; shift 2 ;;
    --rmv-qual) RMV_QUAL="$2"; shift 2 ;;
    --keep-intermediates) KEEP_INTERMEDIATES=true; shift ;;
    --overwrite) OVERWRITE=true; shift ;;
    *) usage ;;
  esac
done

# Validate Mandatory Arguments
if [[ -z "$INPUT_COUNTS" ]]; then echo "Error: -I <hdf5s> is a mandatory flag." >&2; usage; fi
if [[ -z "$OUTPUT_DIR" ]]; then echo "Error: -O <output_dir> is a mandatory flag." >&2; usage; fi
if [[ -z "$REF" ]]; then echo "Error: -R <reference.fa> is a mandatory flag." >&2; usage; fi
# -C is checked further down, once the mode and the input samples are both known, because CASE mode
# can derive it.

# Validate the reference and its companion index files
if [[ ! -f "$REF" ]]; then echo "Error: Reference FASTA not found: ${REF}" >&2; exit 1; fi
if [[ ! -f "${REF}.fai" ]]; then
  echo "Error: Reference index not found: ${REF}.fai. Create it with 'samtools faidx ${REF}'." >&2
  exit 1
fi
REF_DICT="${REF%.*}.dict"
if [[ ! -f "$REF_DICT" ]]; then
  echo "Error: Sequence dictionary not found: ${REF_DICT}. Create it with 'gatk CreateSequenceDictionary -R ${REF}'." >&2
  exit 1
fi

# Validate numeric arguments
for ARG_PAIR in "--min-qual:${MIN_QUAL}" "--rmv-qual:${RMV_QUAL}"; do
  ARG_NAME="${ARG_PAIR%%:*}"
  ARG_VALUE="${ARG_PAIR#*:}"
  if ! [[ "$ARG_VALUE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: ${ARG_NAME} must be a non-negative number (e.g. 100.0)." >&2
    exit 1
  fi
done
for ARG_PAIR in "--low-count-pct:${LOW_COUNT_PCT}" "--scatters:${SCATTERS_REQUESTED}"; do
  ARG_NAME="${ARG_PAIR%%:*}"
  ARG_VALUE="${ARG_PAIR#*:}"
  if ! [[ "$ARG_VALUE" =~ ^[0-9]+$ ]]; then
    echo "Error: ${ARG_NAME} must be a non-negative integer." >&2
    exit 1
  fi
done
if [ "$LOW_COUNT_PCT" -gt 100 ]; then
  echo "Error: --low-count-pct must be a percentage between 0 and 100." >&2
  exit 1
fi
if [ "$SCATTERS_REQUESTED" -lt 1 ]; then
  echo "Error: --scatters must be at least 1." >&2
  exit 1
fi
# The two QUAL thresholds are only sensible if the drop threshold sits below the tag threshold.
if awk -v a="$RMV_QUAL" -v b="$MIN_QUAL" 'BEGIN { exit !(a > b) }'; then
  echo "Error: --rmv-qual (${RMV_QUAL}) must not be higher than --min-qual (${MIN_QUAL})." >&2
  exit 1
fi

# Validate model directory for CASE mode.
CASE_MODE=false
if [[ -n "$MODEL_DIR" ]]; then
  CASE_MODE=true
  MODEL_DIR="${MODEL_DIR%/}"
  if [[ ! -d "$MODEL_DIR" ]]; then echo "Error: Model directory not found: ${MODEL_DIR}" >&2; exit 1; fi
  if [[ ! -d "${MODEL_DIR}/ploidy-model" ]]; then
    echo "Error: No ploidy-model directory inside ${MODEL_DIR}." >&2
    echo "       -M expects the <out_dir>/gcnv/<cohort> directory of an earlier COHORT run." >&2
    exit 1
  fi
  if ! compgen -G "${MODEL_DIR}/gcnvcaller_scatters/scatter_*-model" >/dev/null; then
    echo "Error: No scatter model shards inside ${MODEL_DIR}/gcnvcaller_scatters/." >&2
    echo "       -M expects the <out_dir>/gcnv/<cohort> directory of an earlier COHORT run." >&2
    exit 1
  fi
  if [[ -n "$INTERVAL_FILE" ]]; then
    echo "WARNING: -L is ignored in CASE mode. Intervals come from the model in ${MODEL_DIR}."
  fi
fi

# --- Resolve input read counts ---
if [ "$CASE_MODE" = true ]; then
  if [[ "$INPUT_COUNTS" == *,* ]]; then
    echo "Error: CASE mode calls a single sample, but -I lists several files." >&2
    exit 1
  fi
  COUNT_PATHS=("$INPUT_COUNTS")
elif [ -f "$INPUT_COUNTS" ]; then
  echo "INFO: Reading read count paths from sample sheet ${INPUT_COUNTS}."
  mapfile -t COUNT_PATHS < <(grep -v '^[[:space:]]*$' "$INPUT_COUNTS" | tr -d '\r' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
else
  echo "INFO: Reading read count paths from the command line."
  IFS=',' read -r -a COUNT_PATHS <<< "$INPUT_COUNTS"
fi

COUNT_ARGS=()
for COUNTS in "${COUNT_PATHS[@]}"; do
  if [[ ! -f "$COUNTS" ]]; then echo "Error: Read count file not found: ${COUNTS}" >&2; exit 1; fi
  if [[ "$COUNTS" != *.hdf5 ]]; then
    echo "Error: -I expects the .hdf5 files written by varwolf.sh --counts. Got: ${COUNTS}" >&2
    exit 1
  fi
  COUNT_ARGS+=("-I" "$COUNTS")
done
NUM_SAMPLES=${#COUNT_PATHS[@]}

# --- Sample counts and the name of this run ---

if [ "$CASE_MODE" = false ]; then
  if [[ -z "$COHORT" ]]; then
    echo "Error: -C <cohort_name> is a mandatory flag in COHORT mode. It names the model being built." >&2
    usage
  fi
  if [ "$NUM_SAMPLES" -lt 10 ]; then
    echo "Error: gCNV requires at least 10 samples" >&2
    exit 1
  fi
else
  if [[ -n "$COHORT" ]]; then
    echo "WARNING: -C is ignored in CASE mode. This run is named after the sample."
  fi
  SAMPLE="$(basename "${COUNT_PATHS[0]}" .hdf5)"
fi

# --- Resource sets ---

case "$GENOME_BUILD" in
  hg38)
    ALLOSOMAL=(chrX chrY)
    ;;
  *)
    echo "Error: --build must be 'hg38'. Got '${GENOME_BUILD}'." >&2
    usage
    ;;
esac

# --- Start script ---

if [ "$CASE_MODE" = false ]; then
  GCNV_DIR="${OUTPUT_DIR}/gcnv/${COHORT}"
else
  GCNV_DIR="${OUTPUT_DIR}/gcnv/${SAMPLE}"
fi

# Check earlier run overwriting
if [ "$CASE_MODE" = true ] && [ "$(realpath -m "${GCNV_DIR}")" = "$(realpath -m "${MODEL_DIR}")" ]; then
  echo "Error: This CASE run would work in ${GCNV_DIR}, which is the -M model directory itself." >&2
  echo "       Use a different -O, or a sample whose name differs from the cohort name." >&2
  exit 1
fi
if [ -d "${GCNV_DIR}/ploidy-model" ] || [ -d "${GCNV_DIR}/ploidy-calls" ] \
  || [ -d "${GCNV_DIR}/interval_scatters" ] || [ -d "${GCNV_DIR}/gcnvcaller_scatters" ]; then
  if [ "$OVERWRITE" = false ]; then
    echo "Error: ${GCNV_DIR} already holds results from an earlier run." >&2
    echo "       Re-run with --overwrite to delete them and start over, or choose a different name." >&2
    exit 1
  fi
  echo "WARNING: --overwrite given. Deleting the earlier run in ${GCNV_DIR}..."
  rm -rf "${GCNV_DIR}/ploidy-model" "${GCNV_DIR}/ploidy-calls" \
    "${GCNV_DIR}/interval_scatters" "${GCNV_DIR}/gcnvcaller_scatters"
  rm -f "${GCNV_DIR}/${COHORT}_bins.interval_list" "${GCNV_DIR}/annotated.interval_list" \
    "${GCNV_DIR}/unrestricted.interval_list" "${GCNV_DIR}/filtered.interval_list" \
    "${GCNV_DIR}"/*_denoised_copy_ratios.tsv \
    "${GCNV_DIR}"/*_intervals.cnv.vcf.gz "${GCNV_DIR}"/*_intervals.cnv.vcf.gz.tbi
fi

mkdir -p "${OUTPUT_DIR}/vcfs"
mkdir -p "${GCNV_DIR}"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate menagerie
set -u

# Write PAR and CPP
if [[ -n "$CUSTOM_PAR" ]]; then
  if [[ ! -f "$CUSTOM_PAR" ]]; then echo "Error: PAR BED file not found: ${CUSTOM_PAR}" >&2; exit 1; fi
  PAR_BED="$CUSTOM_PAR"
  echo "INFO: Using custom PAR intervals from ${PAR_BED}."
else
  PAR_BED="${GCNV_DIR}/par.bed"
  printf 'chrX\t10001\t2781479\nchrX\t155701383\t156030895\nchrY\t10001\t2781479\nchrY\t56887903\t57217415\n' \
    > "${PAR_BED}"
  echo "INFO: Wrote ${GENOME_BUILD} PAR intervals to ${PAR_BED}."
fi

if [[ -n "$PLOIDY_PRIORS" ]]; then
  if [[ ! -f "$PLOIDY_PRIORS" ]]; then echo "Error: Ploidy priors file not found: ${PLOIDY_PRIORS}" >&2; exit 1; fi
  echo "INFO: Using custom contig ploidy priors from ${PLOIDY_PRIORS}."
else
  PLOIDY_PRIORS="${GCNV_DIR}/contig_ploidy_priors.tsv"
  {
    printf 'CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3\n'
    for CHR in chr{1..22}; do
      printf '%s\t0.0\t0.01\t0.98\t0.01\n' "$CHR"
    done
    printf 'chrX\t0.01\t0.49\t0.49\t0.01\n'
    printf 'chrY\t0.495\t0.495\t0.01\t0.0\n'
  } > "${PLOIDY_PRIORS}"
  echo "INFO: Wrote ${GENOME_BUILD} contig ploidy priors to ${PLOIDY_PRIORS}."
fi

if [ "$CASE_MODE" = false ]; then
  echo "INFO: Fitting a new gCNV model for cohort '${COHORT}' from ${NUM_SAMPLES} samples."
else
  echo "INFO: Calling CNVs in ${SAMPLE} against the model in ${MODEL_DIR}."
fi

CPUS="${SLURM_CPUS_PER_TASK:-1}"

# --- Interval preparation ---

if [ "$CASE_MODE" = false ]; then
  BINS="${GCNV_DIR}/${COHORT}_bins.interval_list"
  if [[ -n "$INTERVAL_FILE" ]]; then
    echo "INFO: Preprocessing provided intervals..."
    gatk PreprocessIntervals \
      -R "${REF}" \
      -L "${INTERVAL_FILE}" \
      --bin-length 0 \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O "${BINS}"
  else
    echo "INFO: No interval file provided. Creating whole-genome bins..."
    gatk PreprocessIntervals \
      -R "${REF}" \
      --bin-length 1000 \
      --padding 0 \
      -imr OVERLAPPING_ONLY \
      -O "${BINS}"
  fi
  echo "INFO: Created bins!"

  echo "INFO: Annotating intervals..."
  gatk AnnotateIntervals \
    -R "${REF}" \
    -L "${BINS}" \
    -imr OVERLAPPING_ONLY \
    -O "${GCNV_DIR}/annotated.interval_list"
  echo "INFO: Finished annotating intervals!"

  echo "INFO: Filtering intervals..."
  gatk FilterIntervals \
    -L "${BINS}" \
    -XL "${PAR_BED}" \
    --annotated-intervals "${GCNV_DIR}/annotated.interval_list" \
    -imr OVERLAPPING_ONLY \
    "${COUNT_ARGS[@]}" \
    --low-count-filter-percentage-of-samples "${LOW_COUNT_PCT}" \
    -O "${GCNV_DIR}/unrestricted.interval_list"

  echo "INFO: Restricting intervals to the contigs in the ploidy priors table..."
  awk '
    NR == FNR { if (FNR > 1) { keep[$1] = 1 } ; next }
    /^@/ { print ; next }
    ($1 in keep) { print ; next }
    { dropped[$1] = 1 }
    END { for (contig in dropped) { print "INFO: Dropped contig " contig > "/dev/stderr" } }
  ' "${PLOIDY_PRIORS}" "${GCNV_DIR}/unrestricted.interval_list" > "${GCNV_DIR}/filtered.interval_list"

  INTERVAL_COUNT="$(grep -cv '^@' "${GCNV_DIR}/filtered.interval_list" || true)"
  if [ "${INTERVAL_COUNT}" -eq 0 ]; then
    echo "Error: No intervals left after filtering." >&2
    echo "       The contig names in ${PLOIDY_PRIORS} most likely do not match the reference." >&2
    exit 1
  fi
  echo "INFO: Finished filtering intervals! ${INTERVAL_COUNT} intervals remain."

# Scatter handling

  SCATTER_CONTENT=$(( INTERVAL_COUNT / SCATTERS_REQUESTED ))

  echo "INFO: Scattering ${INTERVAL_COUNT} intervals into ${SCATTERS_REQUESTED} shard(s) of ${SCATTER_CONTENT}..."
  gatk IntervalListTools \
    -I "${GCNV_DIR}/filtered.interval_list" \
    -O "${GCNV_DIR}/interval_scatters" \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT "${SCATTER_CONTENT}"
  echo "INFO: Finished scattering intervals!"
fi

SCATTERS=()
if [ "$CASE_MODE" = false ]; then
  for SCATTER_DIR in "${GCNV_DIR}"/interval_scatters/temp_*_of_*; do
    SCATTERS+=("$(basename "${SCATTER_DIR}" | cut -d "_" -f 2)")
  done
else
  for SHARD_DIR in "${MODEL_DIR}"/gcnvcaller_scatters/scatter_*-model; do
    SCATTER="$(basename "${SHARD_DIR}")"
    SCATTER="${SCATTER#scatter_}"
    SCATTERS+=("${SCATTER%-model}")
  done
fi

NUM_SCATTERS=${#SCATTERS[@]}

if [ "$NUM_SCATTERS" -eq 0 ]; then
  echo "Error: No interval shards to call. Nothing to do." >&2
  exit 1
fi

THREADS_PER_JOB=$(( CPUS / NUM_SCATTERS ))
if [ "$THREADS_PER_JOB" -lt 1 ]; then
  THREADS_PER_JOB=1
fi
export OMP_NUM_THREADS="${THREADS_PER_JOB}"
echo "INFO: Running ${NUM_SCATTERS} interval shard(s) with ${THREADS_PER_JOB} thread(s) each."

# --- Contig Ploidy ---
echo "INFO: Determining contig ploidy..."
if [ "$CASE_MODE" = false ]; then
  gatk DetermineGermlineContigPloidy \
    -L "${GCNV_DIR}/filtered.interval_list" \
    -imr OVERLAPPING_ONLY \
    "${COUNT_ARGS[@]}" \
    -O "${GCNV_DIR}/" \
    --output-prefix ploidy \
    --contig-ploidy-priors "${PLOIDY_PRIORS}"
else
  gatk DetermineGermlineContigPloidy \
    --model "${MODEL_DIR}/ploidy-model/" \
    "${COUNT_ARGS[@]}" \
    -O "${GCNV_DIR}/" \
    --output-prefix ploidy
fi
echo "INFO: Finished determining ploidy!"

# --- CNV Calling ---

echo "INFO: Running GermlineCNVCaller per interval shard..."
if [ "$CASE_MODE" = false ]; then
  printf '%s\n' "${SCATTERS[@]}" | xargs -d '\n' -I{} -P "${CPUS}" gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L "${GCNV_DIR}/interval_scatters/temp_{}_of_${NUM_SCATTERS}/scattered.interval_list" \
    --annotated-intervals "${GCNV_DIR}/annotated.interval_list" \
    -imr OVERLAPPING_ONLY \
    "${COUNT_ARGS[@]}" \
    -O "${GCNV_DIR}/gcnvcaller_scatters" \
    --output-prefix "scatter_{}" \
    --contig-ploidy-calls "${GCNV_DIR}/ploidy-calls"
else
  printf '%s\n' "${SCATTERS[@]}" | xargs -d '\n' -I{} -P "${CPUS}" gatk GermlineCNVCaller \
    --run-mode CASE \
    --model "${MODEL_DIR}/gcnvcaller_scatters/scatter_{}-model" \
    "${COUNT_ARGS[@]}" \
    -O "${GCNV_DIR}/gcnvcaller_scatters" \
    --output-prefix "scatter_{}" \
    --contig-ploidy-calls "${GCNV_DIR}/ploidy-calls"
fi
echo "INFO: All interval shards finished!"

# --- Per sample postprocessing ---
export OMP_NUM_THREADS=1
if [ "$CASE_MODE" = false ]; then
  MODEL_ROOT="${GCNV_DIR}"
else
  MODEL_ROOT="${MODEL_DIR}"
fi

MODEL_ARGS=()
CALL_ARGS=()
for SCATTER in "${SCATTERS[@]}"; do
  MODEL_ARGS+=("--model-shard-path" "${MODEL_ROOT}/gcnvcaller_scatters/scatter_${SCATTER}-model")
  CALL_ARGS+=("--calls-shard-path" "${GCNV_DIR}/gcnvcaller_scatters/scatter_${SCATTER}-calls")
done

ALLOSOMAL_ARGS=()
for CONTIG in "${ALLOSOMAL[@]}"; do
  ALLOSOMAL_ARGS+=("--allosomal-contig" "${CONTIG}")
done

SAMPLE_NAMES=()
for i in $(seq 0 $((NUM_SAMPLES - 1))); do
  NAME_FILE="${GCNV_DIR}/gcnvcaller_scatters/scatter_${SCATTERS[0]}-calls/SAMPLE_${i}/sample_name.txt"
  if [[ ! -f "$NAME_FILE" ]]; then
    echo "Error: ${NAME_FILE} is missing." >&2
    exit 1
  fi
  SAMPLE_NAMES+=("$(tr -d '\r\n' < "${NAME_FILE}")")
done

echo "INFO: Postprocessing CNV calls per sample..."
for i in "${!SAMPLE_NAMES[@]}"; do
  SAMPLE="${SAMPLE_NAMES[$i]}"
  printf '%s\n' \
    --sample-index "${i}" \
    --output-genotyped-intervals "${GCNV_DIR}/${SAMPLE}_intervals.cnv.vcf.gz" \
    --output-genotyped-segments "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz" \
    --output-denoised-copy-ratios "${GCNV_DIR}/${SAMPLE}_denoised_copy_ratios.tsv"
done | xargs -d '\n' -n 8 -P "${CPUS}" gatk PostprocessGermlineCNVCalls \
  "${MODEL_ARGS[@]}" \
  "${CALL_ARGS[@]}" \
  --contig-ploidy-calls "${GCNV_DIR}/ploidy-calls/" \
  "${ALLOSOMAL_ARGS[@]}" \
  --sequence-dictionary "${REF_DICT}"
echo "INFO: All samples postprocessed!"

# --- Filtering ---
echo "INFO: Filtering CNV calls..."
printf '%s\n' "${SAMPLE_NAMES[@]}" | xargs -d '\n' -I{} -P "${CPUS}" gatk VariantFiltration \
  -V "${OUTPUT_DIR}/vcfs/{}_raw.cnv.vcf.gz" \
  -filter "QUAL < ${MIN_QUAL}" --filter-name "CNVQUAL" \
  -filter "QUAL < ${RMV_QUAL}" --filter-name "CNVRMV" \
  -O "${OUTPUT_DIR}/vcfs/{}_filtered.cnv.vcf.gz"
echo "INFO: Finished filtering CNV calls!"

# Drop the CNVRMV segments and the reference calls, and write SVTYPE into the INFO field
echo "INFO: Writing final CNV VCFs..."
for SAMPLE in "${SAMPLE_NAMES[@]}"; do
  zgrep -P -v "CNVRMV|N\t\." "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz" \
    | sed 's/\tEND/\tSVTYPE=CNV;END/g' \
    | bgzip -o "${OUTPUT_DIR}/vcfs/${SAMPLE}.cnv.vcf.gz"
  tabix -f "${OUTPUT_DIR}/vcfs/${SAMPLE}.cnv.vcf.gz"
done

# --- Cleanup ---
if [ "$KEEP_INTERMEDIATES" = false ]; then
  for SAMPLE in "${SAMPLE_NAMES[@]}"; do
    rm -f "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz" "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz.tbi" \
      "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz" "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz.tbi"
  done
  rm -f "${GCNV_DIR}/unrestricted.interval_list"
fi

echo "SUCCESS"
