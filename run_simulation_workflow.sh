#!/bin/bash
# Run the MOSAIC simulation workflow (Snakefile_simulation) on SLURM: one job per step and replicate, with the
# threads / memory / time limits defined in the workflow. Run from the repository root inside the
# mic_data_integration environment (or set MOSAIC_ENV to the conda/mamba environment to activate).
#
#   bash run_simulation_workflow.sh                                          # scenarios in run_scenarios
#   bash run_simulation_workflow.sh --config run_scenarios=[zeroinf90] iterations=100
#   bash run_simulation_workflow.sh -n                                       # dry run
# Extra arguments are passed to snakemake. Adjust PARTITION / SBATCH_EXTRA for your cluster.
set -euo pipefail
cd "$(dirname "$0")"
PARTITION=${PARTITION:-scu-cpu}
SBATCH_EXTRA=${SBATCH_EXTRA:-"--exclude=c12,c41"}
MAX_JOBS=${MAX_JOBS:-200}
if [ -n "${MOSAIC_ENV:-}" ]; then source ~/.bashrc; mamba activate "$MOSAIC_ENV"; fi
mkdir -p logs/slurm
snakemake -s Snakefile_simulation --jobs "$MAX_JOBS" --latency-wait 120 --keep-going --rerun-incomplete \
    --cluster "sbatch -p $PARTITION $SBATCH_EXTRA -c {threads} --mem={resources.mem_mb} -t {resources.runtime} -J mosaic-{rule} -o logs/slurm/mosaic-{rule}-%j.out" \
    "$@"
