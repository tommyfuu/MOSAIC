import os
import warnings
from pathlib import Path


configfile: "./config_autism.yml"

# global params
dataset_name = config['dataset_name']
src = config['src']
post_integration_results = config['post_integration_outputs']
used_R_methods = config['used_R_methods']
used_python_methods = config['used_Python_methods']
evaluation_outputs = config['evaluation_outputs']

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

rule integrate:
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

