import os
import sys
import pandas as pd
from tqdm import tqdm

# Add 'scripts' and its subdirectories to Python's module search path to allow importing custom modules. 
sys.path.insert(1, 'scripts/python')

from utils import qvals_bh, enrich_mwu

print('Ok, let\'s go!')

# First, extract annotation information from the metatranscriptome functional annotation
# Load functional annotations
functional_annotations = pd.read_table('data/annotation/functional/130/functional_annotation.emapper.annotations')
# Cut off weird characters from the transcript names
functional_annotations['#query'] = functional_annotations['#query'].str.split(".", n=1, expand=True)[0]
# Rename the query_id column
functional_annotations.rename(columns={'#query': 'transcript_id'}, inplace=True)

transcript_info = functional_annotations[['transcript_id', 'Description']]

transcript_info = transcript_info[transcript_info['Description'] != '-']
transcript_info.columns = ['transcript_id', 'annotation']

genera = ['Phaeocystis']

for genus in tqdm(genera, desc="Processing genus"):
    print(genus)
    # Combine transcript ids, annotation info, and module membership data for a given module
    ## Create the required annotation file
    input_dir = os.path.join("data/analysis/WGCNA_130/transcripts_3", genus)

    # Iterate over all modules in folder
    for filename in tqdm(os.listdir(input_dir), desc="Processing modules..."):
        if filename.endswith("_content.txt"):
            input_file = os.path.join(input_dir, filename)
            module_data = pd.read_table(input_file)
            module = filename.split("_")[0]

            # The data should be ordered by modulemembership
            module_data.sort_values(by=[module], inplace=True, ascending=False)

            # Run MWU
            res = enrich_mwu(module_data, transcript_info, min_size=1)

            # Save results
            res.to_csv(f"{input_dir}/{module}_MWU_results.csv", index=False)

    # Filter and combine the results for significantly higher ranking annotations
    results = pd.DataFrame()

    for filename in tqdm(os.listdir(input_dir), desc="Analyzing modules..."):
        if filename.endswith("_MWU_results.csv"):
            input_file = os.path.join(input_dir, filename)
            module = filename.split("_")[0]
            res = pd.read_csv(input_file)
            res = res[(res['q_val'] < 0.10) & (res['estimate'] > 0.5)]
            # Add module name and reorganize columns
            res['Module'] = module
            res = res[['Module', 'annotation', 'estimate', 'p_val', 'q_val', 'U', 'size']]
            results = pd.concat([results, res])

    results.to_csv(f"{input_dir}/{genus}_MWU_results.csv", index=False)