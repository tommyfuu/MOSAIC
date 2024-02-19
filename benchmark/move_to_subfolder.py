import os
from pathlib import Path
import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--option", type=int, default=1, help='1 for simulation, 2 for method running outputs, 3 for metastasis')
parser.add_argument("-d", "--datatype", default = 'count', help='either count or relab')
parser.add_argument("-r", "--related", default = 'no', help='whether the batch effect is related to library size')
parser.add_argument("-v", "--verbose", default = False, help='Boolean for whether to print incomplete iterations')

args = parser.parse_args()
option = args.option

print(args.option)
print(args.verbose)
or_l = [1, 1.25, 1.5]
cond_effect_val_l = [0, 0.25, 0.5, 0.75, 1]
batch_effect_val_l = [0, 0.25, 0.5, 0.75, 1]

current_num = 1

if option == 1:
    # 1. raw simulation files
    source_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_yesrelation_102023'
    dest_dir = '/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_MIDAS_1000_yesrelation_102023_compartmentalized'

    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    for i in range(51, 1051, 50):
        current_num = i-50
        # make dest dir 
        Path(f'{dest_dir}/range{i-50}_{i-1}').mkdir(parents=True, exist_ok=True)
        print(f'{dest_dir}/range{i-50}_{i-1}')
        while current_num < i:
            print(current_num)
            # copy files
            for odds_ratio in or_l:
                for cond_effect_val in cond_effect_val_l:
                    for batch_effect_val in batch_effect_val_l:
                        if cond_effect_val + batch_effect_val <= 1:
                            shutil.copy2(f'{source_dir}/ibd_150_count_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}.csv', f'{dest_dir}/range{i-50}_{i-1}')
                            shutil.copy2(f'{source_dir}/ibd_150_relab_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}.csv', f'{dest_dir}/range{i-50}_{i-1}')
                            shutil.copy2(f'{source_dir}/ibd_150_meta_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}.csv', f'{dest_dir}/range{i-50}_{i-1}')
            current_num += 1
        print(f'{dest_dir}/range{i-50}_{i-1}')
        # print(f"rm range{i-50}_{i}.tar.lz4")
        print(f"tar cvf - range{i-50}_{i-1} | lz4 > range{i-50}_{i-1}.tar.lz4")
        # print(f'rclone copy range{i-50}_{i-1}.tar.lz4 box_wcmc:mic_bc_benchmark_sim_no_backup/ --progress')

elif option == 2:
    # 2. method output files
    print(args)
    datatype = args.datatype
    related = args.related
    source_dir = f'/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023'
    dest_dir = f'/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized'

    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    for i in range(51, 1051, 50):
        current_num = i-50
        # make dest dir 
        Path(f'{dest_dir}/range{i-50}_{i-1}').mkdir(parents=True, exist_ok=True)
        print(f'{dest_dir}/range{i-50}_{i-1}')
        while current_num < i:
            print(current_num)
            # copy files
            for odds_ratio in or_l:
                for cond_effect_val in cond_effect_val_l:
                    for batch_effect_val in batch_effect_val_l:
                        if cond_effect_val + batch_effect_val <= 1:
                            shutil.copytree(f'{source_dir}/out_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}', f'{dest_dir}/range{i-50}_{i-1}/out_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}')
            current_num += 1
        print(f"tar cvf - /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized/range{i-50}_{i-1} | lz4 > /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized/output_{datatype}_{related}_range{current_num}_{i-1}.tar.lz4")
        
        # print(f"mv /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized/output_count_yes_range{i-50}_{i-1}.tar.lz4 /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized/output_{datatype}_{related}_range{i-50}_{i-1}.tar.lz4")
        print(f'rclone copy /athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_output_{datatype}_{related}relation_102023_compartmentalized/output_{datatype}_{related}_range{i-50}_{i-1}.tar.lz4 box_wcmc:output_{datatype}_{related}_mic_bc_benchmark_sim_backup/ --progress')

elif option == 3:
    # 3. evaluation files    
    datatype = args.datatype
    related = args.related
    source_dir = f'/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_eval_{datatype}_{related}relation_102023'
    dest_dir = f'/athena/linglab/scratch/chf4012/simulation_outputs/simulation_data_eval_{datatype}_{related}relation_102023_compartmentalized'

    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    for i in range(51, 1051, 50):
        current_num = i-50
        # make dest dir 
        Path(f'{dest_dir}/range{i-50}_{i-1}').mkdir(parents=True, exist_ok=True)
        print(f'{dest_dir}/range{i-50}_{i-1}')
        while current_num < i:
            print(current_num)
            # copy files
            for odds_ratio in or_l:
                for cond_effect_val in cond_effect_val_l:
                    for batch_effect_val in batch_effect_val_l:
                        if cond_effect_val + batch_effect_val <= 1:
                            shutil.copytree(f'{source_dir}/out_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}', f'{dest_dir}/range{i-50}_{i-1}/out_{odds_ratio}_{cond_effect_val}_{batch_effect_val}_iter_{current_num}')
            current_num += 1
