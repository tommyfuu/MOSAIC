import os
import warnings
from pathlib import Path


configfile: "./config.yml"

# global params
dataset_name = config['dataset_name']
src = config['src']

rule preprocess:
    '''If MicrobiomeHD or CMD, load and preprocess data; if simulate, simulate data.'''
    output:
        # out_count = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_count_data.csv',
        # out_metadata = f'{src}/cleaned_data/{dataset_name}/{dataset_name}_metadata.csv',
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
        # elif config['data_type'] == 'CMD':
        #     'python scripts/preprocess_CMD.py {src} {dataset_name}'
        # elif config['data_type'] == 'simulate':
        #     'python scripts/simulate_data.py {dataset_name}'

            
    