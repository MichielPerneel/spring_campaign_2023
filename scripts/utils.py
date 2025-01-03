import os
import yaml
import json
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

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
        file1 = os.path.join(file_path, sample + "_R1.fastq.gz")
        file2 = os.path.join(file_path, sample + "_R2.fastq.gz")
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

def process_file(file_path):
    """
    Process a single kallisto output file to extract TPM and count information.
    """
    tsv = pd.read_csv(file_path, sep='\t', usecols=['target_id', 'tpm', 'est_counts'])
    tpm = tsv[['target_id', 'tpm']].rename(columns={'tpm': os.path.basename(os.path.dirname(file_path))})
    counts = tsv[['target_id', 'est_counts']].rename(columns={'est_counts': os.path.basename(os.path.dirname(file_path))})
    return tpm.set_index('target_id'), counts.set_index('target_id')

def merge_kallisto_hpc(input_directory, output_directory):
    """
    Combines reads from all samples quantified against the assembly using kallisto into combined counts and TPM files
    for each station. This version is optimized for HPC systems.

    :param input_directory: Directory containing processed samples by kallisto.
    :param output_directory: Directory where the combined files for all samples will be saved.
    """
    # Find all abundance.tsv files within the input directory
    tsvlist = [os.path.join(dp, f) for dp, dn, filenames in os.walk(input_directory) for f in filenames if f == 'abundance.tsv']
    
    # Extract station numbers from file paths
    station_files = {}
    for filepath in tsvlist:
        # Extract station number (here "51" or "130")
        station = os.path.basename(os.path.dirname(filepath)).split('_')[0]
        if station not in station_files:
            station_files[station] = []
        station_files[station].append(filepath)
    
    # Process files for each station
    for station, files in station_files.items():
        print(f"Processing station {station} with {len(files)} samples.")

        # Use ProcessPoolExecutor to process files in parallel
        with ProcessPoolExecutor(max_workers=16) as executor:
            results = list(executor.map(process_file, files))
        
        # Merge all TPM and counts dataframes for the current station
        all_tpm = pd.concat([result[0] for result in results], axis=1)
        all_counts = pd.concat([result[1] for result in results], axis=1)

        # Write output
        station_output_dir = os.path.join(output_directory, station)
        if not os.path.exists(station_output_dir):
            os.makedirs(station_output_dir)
        all_tpm.reset_index().to_csv(os.path.join(station_output_dir, f"{station}_tpm.csv"), index=False)
        all_counts.reset_index().to_csv(os.path.join(station_output_dir, f"{station}_count.csv"), index=False)

        print(f"Finished processing station {station}. Output saved to {station_output_dir}.")

def calculate_mapping_rate(run_info_path):
    """Calculate the mapping rate from a kallisto run_info.json file."""
    with open(run_info_path, 'r') as f:
        run_info = json.load(f)

    total_reads = run_info['n_processed']
    mapped_reads = run_info['n_pseudoaligned']
    mapping_rate = (mapped_reads / total_reads) * 100

    return mapping_rate

def qvals_bh(pvals: np.array):
    return multipletests(pvals, method="fdr_bh")[1]


def enrich_mwu(ordered_list, annotation: pd.DataFrame,
               min_size=5, alternative="two-sided", print_progress=0):
    """
    Check if some annotation information occurs at the top/bottom of an ordered list of identifiers.

    :param ordered_list: List of KEGG KO, Module, transcript ids, etc. associated with a WGCNA module, ordered by some metric (e.g. moduleMembership).
    :param annotation: columns: "KEGG_ko", "annotation", "KEGG_Module", etc.
    :param min_size: the minimum occurrence of values in annotation category (that are also in ordered list)
        that we need in order to actually run the MWU test.
        I. e. smaller gene sets are skipped.
    :param alternative: passed to `mannwhitneyu`.
    :param print_progress: print a progress message every 'n' tests.
    :return: A data frame with columns:
        - annotation: the annotation category that was tested
        - p_val: the two-sided p-value of the MWU test:
          the more a term is located at the top/bottom, the lower the p-value.
        - q_val: the BH-adjusted p-value.
        - U: the MWU test-statistic.
        - size: the size of the overlap between the number of terms belonging to a annotation category and all terms in the list
        - estimate: the probability that a random term of the given category
        is ranked higher than another random term of the ordered list.
        - direction: can take three values: "top"/"mid"/"bot", if the estimate is bigger/equal/smaller than 0.5.
          For example, "top" means a term is located closer to the top of the given list.
    """
    # Get the indices of all identifiers in the ordered list.
    all_indices = list(range(len(ordered_list)))

    # Group the annotation datas by their categories.
    groups = annotation.groupby("annotation")
    n_groups = len(groups)
    results = []
    
    # Loop over each pathway.
    for i, (annotation, group) in enumerate(groups):
        # Print a progress message if print_progress is set and i is a multiple of print_progress
        if print_progress and i % print_progress == 0:
            print("{:6d} / {:6d} ({:05.2f}%)".format(i, n_groups, i/n_groups * 100))

        # Get the transcript ids of the current pathway.
        group_transcript_ids = group.transcript_id.values
        
        # Find the indices of the transcript ids in the ordered list that are also in the current pathway.
        indices_of_annotation = np.flatnonzero(np.in1d(ordered_list, group_transcript_ids))

        # Only perform the MWU test if the number of overlapping genes is at least min_size
        if len(indices_of_annotation) >= min_size:
            # Perform the Mann-Whitney U test between the indices of all genes in the ordered list and
            # the indices of the genes in the current GO category
            u2, p_val = mannwhitneyu(x=all_indices, y=indices_of_annotation, alternative=alternative)
            results.append({"annotation": annotation, "p_val": p_val, "U": u2, "size": len(indices_of_annotation)})

    # Record the results in a dictionary
    df = pd.DataFrame.from_records(results)
    # Add the corrected p-values
    df.insert(column="q_val", value=qvals_bh(df["p_val"]), loc=2)
    df["estimate"] = df["U"]/(df["size"] * len(all_indices))
    df["direction"] = df["estimate"].apply(lambda e: "mid" if e == 0.5 else ("top" if e > 0.5 else "bot"))
    df["_extremity"] = 0.5 - abs(0.5 - df["estimate"])
    df.sort_values(by=["p_val", "_extremity"], inplace=True)
    df.drop(columns="_extremity", inplace=True)
    return df