#! /bin/bash -l

#SBATCH --partition=scu-cpu
#SBATCH --job-name=olduploads
#SBATCH --time=10:30:00
#SBATCH --mem=7G   # memory requested, units available: K,M,G,T
#SBATCH --cpus-per-task=1
#SBATCH --output uploadrn-%j.out
#SBATCH --error uploadrn-%j.err
#SBATCH --mail-user=chf4012@med.cornell.edu
#SBATCH --mail-type=ALL

source ~/.bashrc

mamba activate bc_benchmark

tar cvf - /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/range901_950 | lz4 > /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/output_relab_no_range901_950.tar.lz4
rclone copy /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/output_relab_no_range901_950.tar.lz4 box_wcmc:output_relab_no_mic_bc_benchmark_sim_backup/ --progress
tar cvf - /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/range951_1000 | lz4 > /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/output_relab_no_range951_1000.tar.lz4
rclone copy /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_relab_norelation_102023_compartmentalized/output_relab_no_range951_1000.tar.lz4 box_wcmc:output_relab_no_mic_bc_benchmark_sim_backup/ --progress

tar cvf - /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_count_yesrelation_102023_compartmentalized/range951_1000 | lz4 > /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_count_yesrelation_102023_compartmentalized/output_count_yes_range951_1000.tar.lz4
rclone copy /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_count_yesrelation_102023_compartmentalized/output_count_yes_range951_1000.tar.lz4 box_wcmc:output_count_yes_mic_bc_benchmark_sim_backup/ --progress