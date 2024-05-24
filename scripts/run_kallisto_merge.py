# run_kallisto_merge.py
import os
from utils import load_config, merge_kallisto_hpc

# Define input and output directories
config = load_config('config.yaml')
input_directory = os.path.join(config['output_dir'], 'quantification')
output_directory = os.path.join(config['output_dir'], 'quantification')

# Call the function
merge_kallisto_hpc(input_directory, output_directory)