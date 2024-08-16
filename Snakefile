configfile: "./config_simulate.yml"
# make sure this is the first line in the Snakefile

import os
import warnings
from pathlib import Path


# global params
dataset_name = config['dataset_name']
src = config['src']
post_integration_results = config['post_integration_outputs']
used_R_methods = config['used_R_methods']
used_python_methods = config['used_Python_methods']
evaluation_outputs = config['evaluation_outputs']
if 'sim_output_root' in config:
    sim_output_root = config['sim_output_root']
    or_l = config['or_l']
    cond_effect_val_l = config['cond_effect_val_l']
    batch_effect_val_l = config['batch_effect_val_l']
    iteration = config['iter']


# rule all
if config['rules_to_run'] == 'all_rw':
    rule all:
        input:
            out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
            out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv',
            out_check = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_complete_confounding.csv',
            last_R_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_R_methods[-1]}.csv',
            last_python_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_python_methods[-1]}.csv',
            out_runtime = f'{post_integration_results}/{dataset_name}/{dataset_name}_runtime.txt',
            out_summary = f'{evaluation_outputs}/{dataset_name}/global_benchmarking_stats_{dataset_name}.csv',
            out_vis = f'{evaluation_outputs}/{dataset_name}/{dataset_name}_multi_PCOA_both_batch.pdf'
elif config['rules_to_run'] == 'all_simulated':
    rule all:
        input:
            out_count = f'{sim_output_root}/simulate/ibd_150_count_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_relab = f'{sim_output_root}/simulate/ibd_150_relab_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_metadata = f'{sim_output_root}/simulate/ibd_150_meta_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_id_batch = f'{sim_output_root}/simulate/ibd_150_id_batch_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
            out_id_cond = f'{sim_output_root}/simulate/ibd_150_id_cond_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
            last_R_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_R_methods[-1]}.csv',
            last_python_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_python_methods[-1]}.csv',
            out_runtime = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_runtime.txt'
elif config['rules_to_run'] == 'preprocess':
    rule all:
        input:
            out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
            out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv',
            out_check = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_complete_confounding.csv'
elif config['rules_to_run'] == 'integrate_rw':
    rule all:
        input:
            out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
            out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv',
            out_check = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_complete_confounding.csv',
            last_R_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_R_methods[-1]}.csv',
            last_python_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_python_methods[-1]}.csv',
            out_runtime = f'{post_integration_results}/{dataset_name}/{dataset_name}_runtime.txt'
elif config['rules_to_run'] == 'integrate_simulated':
    rule all:
        input:
            out_count = f'{sim_output_root}/simulate/ibd_150_count_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_relab = f'{sim_output_root}/simulate/ibd_150_relab_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_metadata = f'{sim_output_root}/simulate/ibd_150_meta_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
            out_id_batch = f'{sim_output_root}/simulate/ibd_150_id_batch_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
            out_id_cond = f'{sim_output_root}/simulate/ibd_150_id_cond_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
            last_R_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_R_methods[-1]}.csv',
            last_python_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_python_methods[-1]}.csv',
            out_runtime = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_runtime.txt'

# rules
rule preprocess:
    '''If MicrobiomeHD or CMD, load and preprocess data; if simulate, simulate data.'''
    output:
        out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
        out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv',
        out_check = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_complete_confounding.csv'
    run:
        if 'microbiomeHD' in config['dataset_name']:
           shell('python3 ./data/step0_microbiomeHD_precleaning.py -s {src} -d {dataset_name}')
           shell('Rscript ./data/step1_clean_and_prune.R microbiomeHD {dataset_name}')
           shell('python3 data/step2_preprocessing_summarystats.py -d {dataset_name}')
           shell('echo {output.out_check}')
        elif 'CMD' in config['dataset_name']:
           shell('Rscript ./data/step1_clean_and_prune.R CMD {dataset_name}')
           shell('python3 data/step2_preprocessing_summarystats.py -d {dataset_name} -b study_name -c disease')
        else:
            shell('python3 ./data/step0_microbiomeHD_precleaning.py -s {src} -d {dataset_name}')
            shell('Rscript ./data/step1_clean_and_prune.R microbiomeHD {dataset_name}')
            shell('python3 data/step2_preprocessing_summarystats.py -d {dataset_name}')
            shell('echo {output.out_check}')

rule simulate:
    '''If simulate, simulate data.'''
    output:
        out_count = f'{sim_output_root}/simulate/ibd_150_count_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_relab = f'{sim_output_root}/simulate/ibd_150_relab_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_metadata = f'{sim_output_root}/simulate/ibd_150_meta_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_id_batch = f'{sim_output_root}/simulate/ibd_150_id_batch_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
        out_id_cond = f'{sim_output_root}/simulate/ibd_150_id_cond_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt'
    run:
        shell('Rscript ./benchmark/generate_data_MIDASim_wrapper.R')

rule integrate_rw:
    '''Integrate data with benchmarking methods.'''
    input:
        out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
        out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_meta_data.csv',
        out_check = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_complete_confounding.csv',
    output:
        last_R_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_R_methods[-1]}.csv',
        last_python_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_python_methods[-1]}.csv',
        out_runtime = f'{post_integration_results}/{dataset_name}/{dataset_name}_runtime.txt'
    run:
        shell('Rscript ./benchmark/methods_benchmarking.R')
        shell('python3 ./benchmark/evaluate.py -o 4')

rule integrate_simulated:
    """Integrate simulated data."""
    input:
        out_count = f'{sim_output_root}/simulate/ibd_150_count_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_relab = f'{sim_output_root}/simulate/ibd_150_relab_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_metadata = f'{sim_output_root}/simulate/ibd_150_meta_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.csv',
        out_id_batch = f'{sim_output_root}/simulate/ibd_150_id_batch_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt',
        out_id_cond = f'{sim_output_root}/simulate/ibd_150_id_cond_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}.txt'
    output:
        last_R_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_R_methods[-1]}.csv',
        last_python_out_integrated = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_{used_python_methods[-1]}.csv',
        out_runtime = f'{sim_output_root}/benchmark/relab_yesrelation/out_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}/ibd_{or_l[-1]}_{cond_effect_val_l[-1]}_{batch_effect_val_l[-1]}_iter_{iteration}_runtime.txt'
    run:
        shell('Rscript ./benchmark/methods_benchmarking_sim_wrapper.R')
        for datatype in ['count', 'relab']:
            for relation in ['yesrelation', 'norelation']:
                shell('python3 ./benchmark/evaluate.py -o 1 -d -i {iteration} -d {datatype} -r {relation}')

            
rule evaluate:
    '''Evaluate the integrated data.'''
    input:
        last_R_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_R_methods[-1]}.csv',
        last_python_out_integrated = f'{post_integration_results}/{dataset_name}/{dataset_name}_{used_python_methods[-1]}.csv',
        out_runtime = f'{post_integration_results}/{dataset_name}/{dataset_name}_runtime.txt',
    output:
        out_summary = f'{evaluation_outputs}/{dataset_name}/global_benchmarking_stats_{dataset_name}.csv',
        out_vis = f'{evaluation_outputs}/{dataset_name}/{dataset_name}_multi_PCOA_both_batch.pdf'
    run:
        shell('python3 ./benchmark/evaluate.py -o 5')

