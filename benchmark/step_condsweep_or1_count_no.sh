#!/bin/bash
# Condition-effect sweep at OR = 1, count data, no batch-library size confounding (Fig 2 panel d for the
# OR 1 / no-confounding main figure). Outputs go next to the OR 1.25 sweep as out_1_<cond>_0.5_iter_<i>.
# Submits generation -> R methods (covar + nocovar) -> Python methods -> evaluation (regular + nocovar RF),
# each step depending on the same iteration of the previous one. Run from benchmark/.
# Usage: bash step_condsweep_or1_count_no.sh [first_iter] [last_iter]
A=${1:-1}; Z=${2:-1000}
ENV='source ~/.bashrc; mamba activate mic_data_integration;'
COMMON="--exclude=c12,c41 -p scu-cpu --qos=normal --array=$A-$Z"
G=$(sbatch --parsable $COMMON -J cs1gen -c 1 --mem=6G -t 2:00:00 -o ../logs/slurm/cs1gen-%A_%a.out \
    --wrap "$ENV Rscript generate_data_MIDASim_condsweep.R \$SLURM_ARRAY_TASK_ID 1 no")
R=$(sbatch --parsable $COMMON -J cs1R -c 10 --mem=20G -t 12:00:00 -o ../logs/slurm/cs1R-%A_%a.out --dependency=aftercorr:$G \
    --wrap "$ENV Rscript methods_benchmarking_sim_condsweep.R \$SLURM_ARRAY_TASK_ID count_no_nocovar_only 1")
P=$(sbatch --parsable $COMMON -J cs1py -c 1 --mem=10G -t 3:00:00 -o ../logs/slurm/cs1py-%A_%a.out --dependency=aftercorr:$G \
    --wrap "$ENV python3 evaluate.py -o 6 -i \$SLURM_ARRAY_TASK_ID -d count --condsweep-only --condsweep-or 1 --condsweep-relation no_relation")
E=$(sbatch --parsable $COMMON -J cs1eval -c 1 --mem=10G -t 4:00:00 -o ../logs/slurm/cs1eval-%A_%a.out --dependency=aftercorr:$R,aftercorr:$P \
    --wrap "$ENV python3 evaluate.py -o 7 -i \$SLURM_ARRAY_TASK_ID -d count --condsweep-only --condsweep-or 1 --condsweep-relation no_relation && python3 evaluate.py -o 7 -i \$SLURM_ARRAY_TASK_ID -d count --condsweep-only --condsweep-or 1 --condsweep-relation no_relation --rerun-rf")
echo "condsweep OR1 count no-relation iters $A-$Z: gen=$G R=$R py=$P eval=$E"
