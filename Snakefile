import os
import warnings
from pathlib import Path


configfile: "./config.yml"

# if config['rules_to_run'] == 'all':

rule preprocess_or_simulate:
    '''If MicrobiomeHD or CMD, load and preprocess data; if simulate, simulate data.'''
    input:
        if config['data_type'] == 'microbiomeHD':
            src = config['src']
        elif config['data_type'] == 'CMD':
            src = ''
        elif config['data_type'] == 'simulate':
            src = ''
        else:
            raise ValueError('data_type must be MicrobiomeHD, CMD, or simulate')
    output:
        if config['data_type'] in ['MicrobiomeHD', 'CMD']:
            dataset_name = config['dataset_name']
    run:
        if config['data_type'] == 'MicrobiomeHD':
           shell('python3 ./data/step0_microbiomeHD_precleaning.py -s {input.src} -d {output.dataset_name}')
           shell('Rscript ./data/step1_clean_and_prune.R microbiomeHD {output.dataset_name}')
           shell('python3 data/step2_preprocessing_summarystats.py -d {output.dataset_name}')
        elif config['data_type'] == 'CMD':
           shell('Rscript ./data/step1_clean_and_prune.R CMD {output.dataset_name}')
           shell('python3 data/step2_preprocessing_summarystats.py -d {output.dataset_name} -b study_name -c disease')
        # elif config['data_type'] == 'CMD':
        #     'python scripts/preprocess_CMD.py {input.src} {output.dataset_name}'
        # elif config['data_type'] == 'simulate':
        #     'python scripts/simulate_data.py {output.dataset_name}'

            
    