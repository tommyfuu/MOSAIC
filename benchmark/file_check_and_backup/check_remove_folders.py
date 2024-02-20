# read all runtime files across 1000 iterations and all conditions and check if they are complete
import argparse
import os
import shutil

or_l = [1, 1.25, 1.5]
cond_effect_val_l = [0, 0.25, 0.5, 0.75, 1]
batch_effect_val_l = [0, 0.25, 0.5, 0.75, 1]

removed_l = []
removed_l.extend([i for i in range(312,316)])
removed_l.extend([i for i in range(327,360)])

parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, default="/athena/linglab/scratch/chf4012/simulation_outputs")
parser.add_argument("-d", "--datatype", default = 'count', help='either count or relab')
parser.add_argument("-r", "--related", default = 'no', help='whether the batch effect is related to library size')
parser.add_argument("-v", "--verbose", default = False, help='Boolean for whether to print incomplete iterations')

args = parser.parse_args()
print(args.verbose)
dir_path = f'{args.dir}/simulation_data_eval_{args.datatype}_{args.related}relation_102023'
for iteration in removed_l:
    current_iteration_good = True
    for or_val in or_l:
        for cond_effect_val in cond_effect_val_l:
            for batch_effect_val in batch_effect_val_l:
                if cond_effect_val+batch_effect_val <= 1:
                    dirname = f'{dir_path}/out_{or_val}_{cond_effect_val}_{batch_effect_val}_iter_{iteration}/'
                    if os.path.exists(dirname):
                        # recursively remove the content of this dir
                        shutil.rmtree(dirname)
                        print(f'removed {dirname}')
                        continue
                    else:
                        print(f'{dirname} does not exist')
    if not current_iteration_good:
        print(f'iteration {iteration} is incomplete')