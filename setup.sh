#!/bin/bash
#
# Creates the conda environment that every script in this repository depends on.
#
# The pipelines call `conda activate menagerie` themselves and expect every external tool to come
# from that environment. Nothing here relies on the cluster's environment modules, so do not
# `module load` anything before or after running this.
#
# Run this on the login node. It only solves and downloads packages, it does not run any analysis.

set -euo pipefail

# --- Environment contents ---
# Versions are pinned so that every user ends up with the same environment. Bump them here.
ENV_NAME="menagerie"

PACKAGES=(
  "python=3.10"          # gatk4 is built against 3.10; plot-vcfstats needs a python too
  "openjdk=17"           # required by GATK 4.6
  "gatk4=4.6.2.0"        # all four pipelines
  "bcftools=1.21"        # sombie.sh, cohorc.sh. 1.9+ required for F_MISSING and annotate --mark-sites
  "samtools=1.21"        # bampire.sh
  "htslib=1.21"          # provides tabix, used by sombie.sh and cohorc.sh
  "bwa-mem2=2.2.1"       # bampire.sh default aligner
  "bwa=0.7.18"           # bampire.sh --legacy-bwa
  "fastp=0.23.4"         # bampire.sh
  "matplotlib-base"      # plot-vcfstats in cohorc.sh renders PNGs with this
)

# --- Usage ---
usage() {
  cat <<EOF
Creates the conda environment required by the menagerie pipelines.

Usage: ./setup.sh [-n <env_name>] [--force] [--verify-only]

Optional flags:
  -n <name>      Name of the environment to create (default: ${ENV_NAME}).
                 The pipeline scripts all activate '${ENV_NAME}', so only change this for testing.
  --force        Delete and recreate the environment if it already exists.
                 Without this, an existing environment is never touched.
  --verify-only  Skip creation and only check that an existing environment has working tools.
  -h, --help     Show this message.

The environment is created under your own conda prefix, so every user runs this once for themselves.
Expect the solve and download to take several minutes.
EOF
  exit 1
}

FORCE=false
VERIFY_ONLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n) ENV_NAME="$2"; shift 2 ;;
    --force) FORCE=true; shift ;;
    --verify-only) VERIFY_ONLY=true; shift ;;
    -h|--help) usage ;;
    *) echo "Error: unknown argument '$1'" >&2; usage ;;
  esac
done

# --- Locate conda ---
if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda is not on your PATH." >&2
  echo "On this cluster it lives at /cm/shared/apps/miniconda3. Add it to your shell with:" >&2
  echo "  /cm/shared/apps/miniconda3/bin/conda init bash && exec bash" >&2
  exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Prefer mamba: the solve for this package set is slow under classic conda.
SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
  echo "INFO: Using mamba to solve the environment."
else
  echo "INFO: mamba not found, falling back to conda. The solve may take a while."
fi

# --- Create the environment ---
if [ "$VERIFY_ONLY" = false ]; then
  if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    if [ "$FORCE" = true ]; then
      echo "INFO: Environment '${ENV_NAME}' exists. Removing it because --force was given..."
      conda env remove -n "${ENV_NAME}" -y
    else
      echo "Error: a conda environment named '${ENV_NAME}' already exists." >&2
      echo "Re-run with --force to delete and recreate it, or with --verify-only to just check it." >&2
      exit 1
    fi
  fi

  echo "INFO: Creating environment '${ENV_NAME}'. This will take several minutes..."
  "${SOLVER}" create -y -n "${ENV_NAME}" \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    "${PACKAGES[@]}"
  echo "INFO: Finished creating '${ENV_NAME}'!"
fi

# --- Verify ---
# Every tool is checked in the environment itself, so a broken install is caught here rather than
# three days into a cohort run.
echo "INFO: Verifying tools in '${ENV_NAME}'..."
conda activate "${ENV_NAME}"

FAILED=()
check() {
  local name="$1"; shift
  if "$@" >/dev/null 2>&1; then
    printf '  OK      %-14s %s\n' "${name}" "$(command -v "${name}")"
  else
    printf '  MISSING %-14s\n' "${name}"
    FAILED+=("${name}")
  fi
}

check gatk          gatk --version
check bcftools      bcftools --version
check samtools      samtools --version
check tabix         tabix --version
check bwa-mem2      bwa-mem2 version
check fastp         fastp --version
check plot-vcfstats plot-vcfstats --help
# bwa has no clean version flag: it exits 1 on --version and prints usage to stderr.
if command -v bwa >/dev/null 2>&1; then
  printf '  OK      %-14s %s\n' "bwa" "$(command -v bwa)"
else
  printf '  MISSING %-14s\n' "bwa"
  FAILED+=("bwa")
fi

# cohorc.sh leans on two bcftools features that older builds do not have, so test them directly
# against a throwaway VCF rather than trusting the version string.
echo "INFO: Checking the bcftools features cohorc.sh depends on..."
TEST_VCF="$(mktemp -t menagerie_check_XXXXXX.vcf)"
trap 'rm -f "${TEST_VCF}"' EXIT
cat > "${TEST_VCF}" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/1	./.
EOF

if bcftools filter -e 'F_MISSING > 0.25' "${TEST_VCF}" >/dev/null 2>&1; then
  echo "  OK      F_MISSING expression"
else
  echo "  BROKEN  F_MISSING expression (bcftools is too old, 1.9+ is required)"
  FAILED+=("bcftools F_MISSING")
fi

if bcftools +fill-tags "${TEST_VCF}" -- -t AC,AN >/dev/null 2>&1; then
  echo "  OK      +fill-tags plugin"
else
  echo "  BROKEN  +fill-tags plugin (BCFTOOLS_PLUGINS may be unset or pointing elsewhere)"
  FAILED+=("bcftools +fill-tags")
fi

echo
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "FAILED: ${#FAILED[@]} check(s) did not pass: ${FAILED[*]}" >&2
  echo "Re-run with --force to rebuild the environment from scratch." >&2
  exit 1
fi

cat <<EOF

SUCCESS

The '${ENV_NAME}' environment is ready. The pipeline scripts activate it themselves, so you do not
need to activate it before submitting a job:

  sbatch -c 64 -o logs/SAMPLE_%x_%A.log bampire.sh -I R1.fastq.gz,R2.fastq.gz -O out -S SAMPLE -R ref.fa

To use the same tools interactively:

  conda activate ${ENV_NAME}
EOF
