# read all runtime files across 1000 iterations and all conditions and check if they are complete
import argparse
import os


or_l = [1, 1.25, 1.5]
cond_effect_val_l = [0, 0.25, 0.5, 0.75, 1]
batch_effect_val_l = [0, 0.25, 0.5, 0.75, 1]


parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, default="/athena/linglab/scratch/chf4012/simulation_outputs")
parser.add_argument("-r", "--related", default = 'no', help='whether the batch effect is related to library size')
parser.add_argument("-v", "--verbose", default = False, help='Boolean for whether to print incomplete iterations')

args = parser.parse_args()
print(args.verbose)

dir_path = f'{args.dir}/simulation_data_MIDAS_1000_{args.related}relation_102023'
for iteration in range(1, 1001):
    if iteration % 50 == 0:
        print(f'iteration {iteration}')

    current_iteration_good = True
    for or_val in or_l:
        for cond_effect_val in cond_effect_val_l:
            for batch_effect_val in batch_effect_val_l:
                if cond_effect_val+batch_effect_val <= 1:
                    filename = f'{dir_path}/ibd_150_relab_{or_val}_{cond_effect_val}_{batch_effect_val}_iter_{iteration}.csv'
                    # print(filename)
                    if not os.path.exists(filename):
                        # print(f'{filename} does not exist')
                        if args.verbose:
                            print(f'{or_val}_{cond_effect_val}_{batch_effect_val}_{iteration} incomplete')
                        current_iteration_good = False
                        continue
    if not current_iteration_good:
        print(f'iteration {iteration} is incomplete')

print("Done!")

# command line eg python3 check_runtime_completeness.py -d relab