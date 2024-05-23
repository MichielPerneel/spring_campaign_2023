import os
import yaml
import pandas as pd

def load_config(config_file):
    with open(config_file) as f:
        return yaml.load(f, Loader=yaml.FullLoader)
            
def extract_samplenames(excel_path):
    """
    Extract sample names from the excel file.

    :param excel_path: Path to the excel file
    :return: List of sample names
    """
    df = pd.read_excel(excel_path)
    return df["sample_name"].tolist()

def check_raw_files(sample_names, file_path):
    """
    Check if all paired-end files exist for every sample.

    :param sample_names: List of sample names
    :param file_path: Path to the directory containing the raw files
    :return: None
    """
    for sample in sample_names:
        file1 = os.path.join(file_path, sample + "_1.fastq.gz")
        file2 = os.path.join(file_path, sample + "_2.fastq.gz")
        if not os.path.isfile(file1):
            raise Exception("File {} does not exist.".format(file1))
        if not os.path.isfile(file2):
            raise Exception("File {} does not exist.".format(file2))

def check_config(config):
    """ This function checks if the config file is valid. 
    
    :param config: Config file
    :return: None
    """
    required_keys = ["home_dir", "samples", "raw_reads",
                     "output_dir", "scratch_dir", "adapters",
                     "conda_environments", "ERCC_folder", "spikefile"]
    for key in required_keys:
        # Check if the config file contains all terms
        if key not in config:
            raise Exception("Key {} is missing in the config file.".format(key))
         # Check if the path exists
        if not os.path.exists(config[key]):
            raise Exception("The path for key {} does not exist: {}".format(key, config[key]))
        
    samplenames = extract_samplenames(config["samples"])
    check_raw_files(samplenames, config["raw_reads"])

# Add the transformed merge_kallisto function
def merge_kallisto(input_directory, output_directory):
    """
    Combines reads from all samples quantified against the assembly using kallisto into a single counts and TPM file.

    :param input_directory: Directory containing processed samples by kallisto.
    :param output_directory: Directory where the combined files for all samples will be saved.
    """
    tsvlist = os.listdir(input_directory)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    out_tpm = os.path.join(output_directory, "tpm.csv")
    out_count = os.path.join(output_directory, "count.csv")

    # Initialize dataframes for TPM and counts
    data_tpm, data_counts = pd.DataFrame(), pd.DataFrame()

    # Iterate over sample output folders
    for i, folder in enumerate(tsvlist):
        tsv_path = os.path.join(input_directory, folder, 'abundance.tsv')
        tsv = pd.read_csv(tsv_path, sep='\t')

        # Process TPM
        tpm = tsv[['target_id', 'tpm']].rename(columns={'tpm': folder})
        if i == 0:
            data_tpm = tpm.set_index('target_id')
        else:
            data_tpm = data_tpm.join(tpm.set_index('target_id'), how='inner')

        # Process counts
        counts = tsv[['target_id', 'est_counts']].rename(columns={'est_counts': folder})
        if i == 0:
            data_counts = counts.set_index('target_id')
        else:
            data_counts = data_counts.join(counts.set_index('target_id'), how='inner')

    # Write output files
    data_tpm.reset_index().to_csv(out_tpm, index=False)
    data_counts.reset_index().to_csv(out_count, index=False)
