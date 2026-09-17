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
  Provide the <out_dir>/cnv/<cohort> directory that the earlier run produced.

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
  -M <dir>       The <out_dir>/cnv/<cohort> directory of an earlier COHORT run for calling in CASE mode.
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
                         (default 10). The intervals left after filtering are divided into this many
                         equal parts, whatever the panel size, and the shards run concurrently.
                         Use --scatters 1 to skip scattering and call every interval in one job.
  --ploidy-priors <file> Contig ploidy priors table for DetermineGermlineContigPloidy. Defaults to a
                         standard diploid autosome table for the selected build, written into the
                         working directory. Supply your own for non-chr contig names, or to call
                         contigs this script skips.
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
  <out_dir>/vcfs/<sample>_intervals.cnv.vcf.gz  Per interval genotypes, kept for inspection.
  <out_dir>/cnv/<cohort>/                       Ploidy and gCNV models, scattered calls and denoised
                                                copy ratios. Never deleted unless --overwrite is given,
                                                so it can be reused as the -M directory for later CASE runs.

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

# Validate CASE mode inputs. The layout checked for here is exactly the one a COHORT run leaves
# behind, so a model directory is either a complete previous run or a mistake.
CASE_MODE=false
if [[ -n "$MODEL_DIR" ]]; then
  CASE_MODE=true
  MODEL_DIR="${MODEL_DIR%/}"
  if [[ ! -d "$MODEL_DIR" ]]; then echo "Error: Model directory not found: ${MODEL_DIR}" >&2; exit 1; fi
  if [[ ! -d "${MODEL_DIR}/ploidy-model" ]]; then
    echo "Error: No ploidy-model directory inside ${MODEL_DIR}." >&2
    echo "       -M expects the <out_dir>/cnv/<cohort> directory of an earlier COHORT run." >&2
    exit 1
  fi
  if ! compgen -G "${MODEL_DIR}/gcnvcaller_scatters/scatter_*-model" >/dev/null; then
    echo "Error: No scatter model shards inside ${MODEL_DIR}/gcnvcaller_scatters/." >&2
    echo "       -M expects the <out_dir>/cnv/<cohort> directory of an earlier COHORT run." >&2
    exit 1
  fi
  if [[ -n "$INTERVAL_FILE" ]]; then
    echo "WARNING: -L is ignored in CASE mode. Intervals come from the model in ${MODEL_DIR}."
  fi
fi

# --- Resolve input read counts ---
# What -I means depends on the mode. A CASE run calls one sample, so -I is that one .hdf5 file and
# there is nothing to resolve. A COHORT run takes many, either as a sample sheet or as a comma
# separated list on the command line. Splitting it this way also avoids a trap: a lone .hdf5 path is
# an existing file, so a plain existence test would read the counts themselves as a sample sheet.
# The inputs are resolved before the run is named, because in CASE mode the name comes from the
# sample itself.
if [ "$CASE_MODE" = true ]; then
  if [[ "$INPUT_COUNTS" == *,* ]]; then
    echo "Error: CASE mode calls a single sample, but -I lists several files." >&2
    echo "       Submit one job per sample, all against the same -M." >&2
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
# The two modes want opposite things from their inputs, and that settles what -C has to mean.
#
# COHORT mode fits a model, which needs a crowd, and -C names the model that comes out.
#
# CASE mode calls one sample against a model somebody else already built and named. There is no
# cohort to name, so the run is named after the sample, the way varwolf.sh derives -S from the BAM,
# and -C is not used at all. Several samples in CASE mode are several independent sbatch jobs
# sharing one model, which is also how they run in parallel across nodes rather than in one job.
if [ "$CASE_MODE" = false ]; then
  if [[ -z "$COHORT" ]]; then
    echo "Error: -C <cohort_name> is a mandatory flag in COHORT mode. It names the model being built." >&2
    usage
  fi
  if [ "$NUM_SAMPLES" -lt 10 ]; then
    echo "Error: Only ${NUM_SAMPLES} read count file(s) provided. Fitting a gCNV model needs a" >&2
    echo "       reasonably large cohort, so this script requires at least 10 samples in COHORT mode." >&2
    echo "       To call a small number of samples, build a model from a larger batch first and then" >&2
    echo "       rerun this script with -M." >&2
    exit 1
  fi
  if [ "$NUM_SAMPLES" -lt 30 ]; then
    echo "WARNING: Only ${NUM_SAMPLES} samples. GATK recommends at least 30 to fit a gCNV model."
    echo "WARNING: Expect noisy calls, especially on the sex chromosomes."
  fi
else
  if [[ -n "$COHORT" ]]; then
    echo "WARNING: -C is ignored in CASE mode. This run is named after the sample."
  fi
  SAMPLE="$(basename "${COUNT_PATHS[0]}" .hdf5)"
fi

# --- Resource sets ---
# Each build carries its pseudoautosomal regions and the contigs that get called. Both are small,
# fixed coordinate sets, so they are written out from here rather than pointed at files on this
# cluster that would have to be copied along with the scripts.
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
# Create output directories. A COHORT run's cnv directory is named after the cohort and holds the
# model, which is never auto-deleted, the same way cohorc.sh keeps its GenomicsDB workspaces, so
# later CASE runs can point -M at it. A CASE run is named after its one sample instead, so that
# every sample called against a shared model keeps its own working directory and the jobs can run
# side by side without treading on each other.
if [ "$CASE_MODE" = false ]; then
  CNV_DIR="${OUTPUT_DIR}/cnv/${COHORT}"
else
  CNV_DIR="${OUTPUT_DIR}/cnv/${SAMPLE}"
fi

# --- Earlier runs ---
# IntervalListTools and GermlineCNVCaller never clear their output directories, so shards left behind
# by an earlier run with the same name get picked up alongside the new ones. With a different shard
# count that means duplicate shard numbers writing into the same directories at once, and in CASE
# mode a stale model shard. So a directory that already holds results is refused unless --overwrite
# says to clear it first.
if [ "$CASE_MODE" = true ] && [ "$(realpath -m "${CNV_DIR}")" = "$(realpath -m "${MODEL_DIR}")" ]; then
  echo "Error: This CASE run would work in ${CNV_DIR}, which is the -M model directory itself." >&2
  echo "       Use a different -O, or a sample whose name differs from the cohort name." >&2
  exit 1
fi
if [ -d "${CNV_DIR}/ploidy-model" ] || [ -d "${CNV_DIR}/ploidy-calls" ] \
  || [ -d "${CNV_DIR}/interval_scatters" ] || [ -d "${CNV_DIR}/gcnvcaller_scatters" ]; then
  if [ "$OVERWRITE" = false ]; then
    echo "Error: ${CNV_DIR} already holds results from an earlier run." >&2
    echo "       Re-run with --overwrite to delete them and start over, or choose a different name." >&2
    exit 1
  fi
  echo "WARNING: --overwrite given. Deleting the earlier run in ${CNV_DIR}..."
  if [ "$CASE_MODE" = false ]; then
    echo "WARNING: Any CASE runs called against the old '${COHORT}' model will no longer match it."
  fi
  rm -rf "${CNV_DIR}/ploidy-model" "${CNV_DIR}/ploidy-calls" \
    "${CNV_DIR}/interval_scatters" "${CNV_DIR}/gcnvcaller_scatters" "${CNV_DIR}/jobs"
  rm -f "${CNV_DIR}/${COHORT}_bins.interval_list" "${CNV_DIR}/annotated.interval_list" \
    "${CNV_DIR}/unrestricted.interval_list" "${CNV_DIR}/filtered.interval_list" \
    "${CNV_DIR}"/*_denoised_copy_ratios.tsv
fi

mkdir -p "${OUTPUT_DIR}/vcfs"
mkdir -p "${CNV_DIR}"

# Initiate conda environment
set +u
eval "$(conda shell.bash hook)"
conda activate menagerie
set -u

# --- Write the build resources ---
# Pseudoautosomal regions. gCNV models chrX and chrY as haploid or diploid per sample, which the PAR
# breaks, so these intervals are excluded from the model entirely.
if [[ -n "$CUSTOM_PAR" ]]; then
  if [[ ! -f "$CUSTOM_PAR" ]]; then echo "Error: PAR BED file not found: ${CUSTOM_PAR}" >&2; exit 1; fi
  PAR_BED="$CUSTOM_PAR"
  echo "INFO: Using custom PAR intervals from ${PAR_BED}."
else
  PAR_BED="${CNV_DIR}/par.bed"
  printf 'chrX\t10001\t2781479\nchrX\t155701383\t156030895\nchrY\t10001\t2781479\nchrY\t56887903\t57217415\n' \
    > "${PAR_BED}"
  echo "INFO: Wrote ${GENOME_BUILD} PAR intervals to ${PAR_BED}."
fi

# Contig ploidy priors. Autosomes are diploid with a little room for whole chromosome events, chrX is
# split between one and two copies, and chrY is either absent or single copy. Contigs missing from
# this table are not called at all, which is what keeps chrM and the alt contigs out of the model.
if [[ -n "$PLOIDY_PRIORS" ]]; then
  if [[ ! -f "$PLOIDY_PRIORS" ]]; then echo "Error: Ploidy priors file not found: ${PLOIDY_PRIORS}" >&2; exit 1; fi
  echo "INFO: Using custom contig ploidy priors from ${PLOIDY_PRIORS}."
else
  PLOIDY_PRIORS="${CNV_DIR}/contig_ploidy_priors.tsv"
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

# --- Interval preparation (COHORT mode only) ---
# In CASE mode all of this is inherited from the model, which is the whole point of CASE mode: the
# new samples have to be called over exactly the intervals the model was fitted on.
if [ "$CASE_MODE" = false ]; then
  # Reproduce the bins varwolf.sh --counts made, so that the interval list lines up with the HDF5s.
  # These two commands must stay identical to the ones in varwolf.sh.
  BINS="${CNV_DIR}/${COHORT}_bins.interval_list"
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

  # GC content per interval. FilterIntervals uses it to drop intervals that gCNV cannot model, and
  # GermlineCNVCaller uses it to correct the coverage it sees.
  echo "INFO: Annotating intervals..."
  gatk AnnotateIntervals \
    -R "${REF}" \
    -L "${BINS}" \
    -imr OVERLAPPING_ONLY \
    -O "${CNV_DIR}/annotated.interval_list"
  echo "INFO: Finished annotating intervals!"

  echo "INFO: Filtering intervals..."
  gatk FilterIntervals \
    -L "${BINS}" \
    -XL "${PAR_BED}" \
    --annotated-intervals "${CNV_DIR}/annotated.interval_list" \
    -imr OVERLAPPING_ONLY \
    "${COUNT_ARGS[@]}" \
    --low-count-filter-percentage-of-samples "${LOW_COUNT_PCT}" \
    -O "${CNV_DIR}/unrestricted.interval_list"

  # DetermineGermlineContigPloidy refuses to run if the intervals cover a contig the priors table
  # does not mention, so the priors table decides what gets called. Dropping the rest here turns
  # what would be a crash several minutes in into a line of log output.
  echo "INFO: Restricting intervals to the contigs in the ploidy priors table..."
  awk '
    NR == FNR { if (FNR > 1) { keep[$1] = 1 } ; next }
    /^@/ { print ; next }
    ($1 in keep) { print ; next }
    { dropped[$1] = 1 }
    END { for (contig in dropped) { print "INFO: Dropped contig " contig > "/dev/stderr" } }
  ' "${PLOIDY_PRIORS}" "${CNV_DIR}/unrestricted.interval_list" > "${CNV_DIR}/filtered.interval_list"

  # --- Scatter ---
  # SCATTER_CONTENT is an interval count per shard, but the number of shards is what actually
  # matters, since it sets how many GermlineCNVCaller processes run at once and how many cores each
  # of them gets. So the requested shard count is what the flag takes, and the per shard interval
  # count is derived from it here, rounding up so the division never leaves a stray remainder shard.
  INTERVAL_COUNT="$(grep -cv '^@' "${CNV_DIR}/filtered.interval_list" || true)"
  if [ "${INTERVAL_COUNT}" -eq 0 ]; then
    echo "Error: No intervals left after filtering." >&2
    echo "       The contig names in ${PLOIDY_PRIORS} most likely do not match the reference." >&2
    exit 1
  fi
  echo "INFO: Finished filtering intervals! ${INTERVAL_COUNT} intervals remain."

  SCATTER_TARGET="${SCATTERS_REQUESTED}"
  if [ "${SCATTER_TARGET}" -gt "${INTERVAL_COUNT}" ]; then
    echo "WARNING: Asked for ${SCATTER_TARGET} shards but there are only ${INTERVAL_COUNT} intervals."
    SCATTER_TARGET="${INTERVAL_COUNT}"
  fi
  SCATTER_CONTENT=$(( (INTERVAL_COUNT + SCATTER_TARGET - 1) / SCATTER_TARGET ))

  echo "INFO: Scattering ${INTERVAL_COUNT} intervals into ${SCATTER_TARGET} shard(s) of ${SCATTER_CONTENT}..."
  gatk IntervalListTools \
    -I "${CNV_DIR}/filtered.interval_list" \
    -O "${CNV_DIR}/interval_scatters" \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT "${SCATTER_CONTENT}"
  echo "INFO: Finished scattering intervals!"
fi

# --- Collect the shards to run ---
# The shards are read off disk rather than counted, so that the numbering always agrees with whatever
# IntervalListTools actually produced, or with whatever the model directory actually contains.
SCATTERS=()
SCATTER_LISTS=()
if [ "$CASE_MODE" = false ]; then
  for SCATTER_DIR in "${CNV_DIR}"/interval_scatters/temp_*_of_*; do
    SCATTERS+=("$(basename "${SCATTER_DIR}" | cut -d "_" -f 2)")
    SCATTER_LISTS+=("${SCATTER_DIR}/scattered.interval_list")
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

# Every shard is a separate python process, and gcnvkernel's linear algebra grabs every core on the
# node unless it is told otherwise. The cores are divided between the shards instead, so that the
# shards do not spend the run fighting each other for the same CPUs.
THREADS_PER_JOB=$(( CPUS / NUM_SCATTERS ))
if [ "$THREADS_PER_JOB" -lt 1 ]; then
  THREADS_PER_JOB=1
fi
export OMP_NUM_THREADS="${THREADS_PER_JOB}"
export MKL_NUM_THREADS="${THREADS_PER_JOB}"
export OPENBLAS_NUM_THREADS="${THREADS_PER_JOB}"
# Never more shards at once than there are cores, in case --scatters is set higher than the thread count.
SHARD_JOBS="${NUM_SCATTERS}"
if [ "$SHARD_JOBS" -gt "$CPUS" ]; then
  SHARD_JOBS="${CPUS}"
fi
echo "INFO: Running ${NUM_SCATTERS} interval shard(s), ${SHARD_JOBS} at a time, with ${THREADS_PER_JOB} thread(s) each."

# --- Parallel jobs ---
# Every parallel step below writes one small bash script per job into JOB_DIR and hands the list to
# xargs, which runs at most -P of them at once. xargs exits non-zero if any single job failed, so a
# dead shard or sample can never be mistaken for a finished run.
# The jobs are written to files rather than passed to xargs as command strings because each command
# carries an -I argument per sample, and xargs refuses command lines over 128 KB, which a large cohort
# would reach. The files also leave the exact failing command on disk to inspect or rerun by hand.
JOB_DIR="${CNV_DIR}/jobs"
rm -rf "${JOB_DIR}"
mkdir -p "${JOB_DIR}"

# --- Contig Ploidy ---
# Ploidy has to be settled before the CNV calling, because gCNV calls copy number relative to the
# ploidy of the contig in that particular sample.
echo "INFO: Determining contig ploidy..."
if [ "$CASE_MODE" = false ]; then
  gatk DetermineGermlineContigPloidy \
    -L "${CNV_DIR}/filtered.interval_list" \
    -imr OVERLAPPING_ONLY \
    "${COUNT_ARGS[@]}" \
    -O "${CNV_DIR}/" \
    --output-prefix ploidy \
    --contig-ploidy-priors "${PLOIDY_PRIORS}"
else
  gatk DetermineGermlineContigPloidy \
    --model "${MODEL_DIR}/ploidy-model/" \
    "${COUNT_ARGS[@]}" \
    -O "${CNV_DIR}/" \
    --output-prefix ploidy
fi
echo "INFO: Finished determining ploidy!"

# --- CNV Calling ---
# printf %q quotes every argument, so paths with spaces or odd characters survive being written into
# the job script and read back by bash.
echo "INFO: Running GermlineCNVCaller per interval shard..."
JOB_FILES=()
for i in "${!SCATTERS[@]}"; do
  SCATTER="${SCATTERS[$i]}"
  JOB_FILE="${JOB_DIR}/gcnvcaller_scatter_${SCATTER}.sh"
  {
    echo "set -euo pipefail"
    if [ "$CASE_MODE" = false ]; then
      printf '%q ' gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L "${SCATTER_LISTS[$i]}" \
        --annotated-intervals "${CNV_DIR}/annotated.interval_list" \
        -imr OVERLAPPING_ONLY \
        "${COUNT_ARGS[@]}" \
        -O "${CNV_DIR}/gcnvcaller_scatters" \
        --output-prefix "scatter_${SCATTER}" \
        --contig-ploidy-calls "${CNV_DIR}/ploidy-calls"
    else
      printf '%q ' gatk GermlineCNVCaller \
        --run-mode CASE \
        --model "${MODEL_DIR}/gcnvcaller_scatters/scatter_${SCATTER}-model" \
        "${COUNT_ARGS[@]}" \
        -O "${CNV_DIR}/gcnvcaller_scatters" \
        --output-prefix "scatter_${SCATTER}" \
        --contig-ploidy-calls "${CNV_DIR}/ploidy-calls"
    fi
    echo
    printf 'echo %q\n' "INFO: Finished GermlineCNVCaller for shard ${SCATTER}!"
  } > "${JOB_FILE}"
  JOB_FILES+=("${JOB_FILE}")
done
if ! printf '%s\n' "${JOB_FILES[@]}" | xargs -d '\n' -n 1 -P "${SHARD_JOBS}" bash; then
  echo "Error: At least one GermlineCNVCaller shard failed. The GATK error is further up in this log." >&2
  echo "       The command for every shard is in ${JOB_DIR}/." >&2
  exit 1
fi
echo "INFO: All interval shards finished!"

# From here on every job handles a single sample and runs single threaded, so up to one job per core
# runs at once, instead of every sample in the cohort starting together.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
SAMPLE_JOBS="${CPUS}"

# --- Per sample postprocessing ---
# PostprocessGermlineCNVCalls stitches the shards back together one sample at a time. The model
# shards come from wherever the model lives, the call shards are always the ones just produced.
if [ "$CASE_MODE" = false ]; then
  MODEL_ROOT="${CNV_DIR}"
else
  MODEL_ROOT="${MODEL_DIR}"
fi

MODEL_ARGS=()
CALL_ARGS=()
for SCATTER in "${SCATTERS[@]}"; do
  MODEL_ARGS+=("--model-shard-path" "${MODEL_ROOT}/gcnvcaller_scatters/scatter_${SCATTER}-model")
  CALL_ARGS+=("--calls-shard-path" "${CNV_DIR}/gcnvcaller_scatters/scatter_${SCATTER}-calls")
done

ALLOSOMAL_ARGS=()
for CONTIG in "${ALLOSOMAL[@]}"; do
  ALLOSOMAL_ARGS+=("--allosomal-contig" "${CONTIG}")
done

# The sample order inside the call shards is the order the HDF5s were given in, but the names are
# read back from the shard rather than derived from the file names, so that the VCFs are named after
# the sample the counts actually belong to.
SAMPLE_NAMES=()
for i in $(seq 0 $((NUM_SAMPLES - 1))); do
  NAME_FILE="${CNV_DIR}/gcnvcaller_scatters/scatter_${SCATTERS[0]}-calls/SAMPLE_${i}/sample_name.txt"
  if [[ ! -f "$NAME_FILE" ]]; then
    echo "Error: Expected ${NUM_SAMPLES} samples in the call shards, but ${NAME_FILE} is missing." >&2
    exit 1
  fi
  SAMPLE_NAMES+=("$(tr -d '\r\n' < "${NAME_FILE}")")
done

echo "INFO: Postprocessing CNV calls per sample, ${SAMPLE_JOBS} at a time..."
JOB_FILES=()
for i in $(seq 0 $((NUM_SAMPLES - 1))); do
  SAMPLE="${SAMPLE_NAMES[$i]}"
  JOB_FILE="${JOB_DIR}/postprocess_${SAMPLE}.sh"
  {
    echo "set -euo pipefail"
    printf '%q ' gatk PostprocessGermlineCNVCalls \
      "${MODEL_ARGS[@]}" \
      "${CALL_ARGS[@]}" \
      --sample-index "${i}" \
      --output-genotyped-intervals "${OUTPUT_DIR}/vcfs/${SAMPLE}_intervals.cnv.vcf.gz" \
      --output-genotyped-segments "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz" \
      --output-denoised-copy-ratios "${CNV_DIR}/${SAMPLE}_denoised_copy_ratios.tsv" \
      --contig-ploidy-calls "${CNV_DIR}/ploidy-calls/" \
      "${ALLOSOMAL_ARGS[@]}" \
      --sequence-dictionary "${REF_DICT}"
    echo
    printf 'echo %q\n' "INFO: Finished postprocessing ${SAMPLE}!"
  } > "${JOB_FILE}"
  JOB_FILES+=("${JOB_FILE}")
done
if ! printf '%s\n' "${JOB_FILES[@]}" | xargs -d '\n' -n 1 -P "${SAMPLE_JOBS}" bash; then
  echo "Error: PostprocessGermlineCNVCalls failed for at least one sample. See further up in this log." >&2
  echo "       The command for every sample is in ${JOB_DIR}/." >&2
  exit 1
fi
echo "INFO: All samples postprocessed!"

# --- Filtering ---
# Both thresholds are recorded in the FILTER column first and only then acted on, so that the tagged
# VCF is a complete record of what was thrown away.
echo "INFO: Filtering CNV calls..."
JOB_FILES=()
for SAMPLE in "${SAMPLE_NAMES[@]}"; do
  JOB_FILE="${JOB_DIR}/filter_${SAMPLE}.sh"
  {
    echo "set -euo pipefail"
    printf '%q ' gatk VariantFiltration \
      -V "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz" \
      -filter "QUAL < ${MIN_QUAL}" --filter-name "CNVQUAL" \
      -filter "QUAL < ${RMV_QUAL}" --filter-name "CNVRMV" \
      -O "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz"
    echo
  } > "${JOB_FILE}"
  JOB_FILES+=("${JOB_FILE}")
done
if ! printf '%s\n' "${JOB_FILES[@]}" | xargs -d '\n' -n 1 -P "${SAMPLE_JOBS}" bash; then
  echo "Error: VariantFiltration failed for at least one sample. See further up in this log." >&2
  echo "       The command for every sample is in ${JOB_DIR}/." >&2
  exit 1
fi
echo "INFO: Finished filtering CNV calls!"

# Drop the CNVRMV segments and the reference calls, and write SVTYPE into the INFO field, which
# PostprocessGermlineCNVCalls declares in the header but leaves out of the records.
echo "INFO: Writing final CNV VCFs..."
JOB_FILES=()
for SAMPLE in "${SAMPLE_NAMES[@]}"; do
  JOB_FILE="${JOB_DIR}/final_vcf_${SAMPLE}.sh"
  {
    echo "set -euo pipefail"
    printf '%q ' zgrep -P -v 'CNVRMV|N\t\.' "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz"
    printf '| '
    printf '%q ' sed 's/\tEND/\tSVTYPE=CNV;END/g'
    printf '| '
    printf '%q ' bgzip -o "${OUTPUT_DIR}/vcfs/${SAMPLE}.cnv.vcf.gz"
    echo
    printf '%q ' tabix -f "${OUTPUT_DIR}/vcfs/${SAMPLE}.cnv.vcf.gz"
    echo
  } > "${JOB_FILE}"
  JOB_FILES+=("${JOB_FILE}")
done
if ! printf '%s\n' "${JOB_FILES[@]}" | xargs -d '\n' -n 1 -P "${SAMPLE_JOBS}" bash; then
  echo "Error: Writing the final CNV VCF failed for at least one sample. See further up in this log." >&2
  echo "       The command for every sample is in ${JOB_DIR}/." >&2
  exit 1
fi

# --- Cleanup ---
# Only the intermediates this script created per sample are removed, and by name rather than by
# glob, so nothing another pipeline left in the vcfs directory is ever caught. The per interval
# genotypes stay, and so does everything under the cnv directory, which is the reusable model.
if [ "$KEEP_INTERMEDIATES" = false ]; then
  for SAMPLE in "${SAMPLE_NAMES[@]}"; do
    rm -f "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz" "${OUTPUT_DIR}/vcfs/${SAMPLE}_raw.cnv.vcf.gz.tbi" \
      "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz" "${OUTPUT_DIR}/vcfs/${SAMPLE}_filtered.cnv.vcf.gz.tbi"
  done
  rm -f "${CNV_DIR}/unrestricted.interval_list"
  rm -rf "${JOB_DIR}"
fi

echo "SUCCESS"
