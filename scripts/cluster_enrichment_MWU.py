import os
import sys
import pandas as pd
from tqdm import tqdm

# Add 'scripts' and its subdirectories to Python's module search path to allow importing custom modules. 
sys.path.insert(1, 'scripts/')

from utils import qvals_bh, enrich_mwu

# First, extract annotation information from the metatranscriptome functional annotation
# Load functional annotations
functional_annotations = pd.read_table('data/annotation/functional_eggnog/130/functional_annotation.emapper.annotations')
# Cut off weird characters from the transcript names
functional_annotations['#query'] = functional_annotations['#query'].str.split(".", n=1, expand=True)[0]
# Rename the query_id column
functional_annotations.rename(columns={'#query': 'transcript_id'}, inplace=True)

transcript_info = functional_annotations[['transcript_id', 'KEGG_ko']]

transcript_info = transcript_info[transcript_info['KEGG_ko'] != '-']
transcript_info.columns = ['transcript_id', 'annotation']

# Combine transcript ids, annotation info, and cluster membership data
clusters = pd.read_csv('data/analysis/identified_clusters_transcripts_correlation.csv')

# Extract the unique cluster names
cluster_names = clusters['Cluster'].unique()

results = pd.DataFrame()

# Iterate over all clusters
for cluster in tqdm(cluster_names, desc="Processing clusters..."):
    # Print cluster name and amount of transcripts in the cluster
    print(f'Processing cluster {cluster} with {clusters[clusters["Cluster"] == cluster].shape[0]} transcripts')
    
    data = clusters[clusters['Cluster'] == cluster].sort_values(by='correlation', ascending=False)
    
    # Add transcript annotation information
    data = pd.merge(data, transcript_info, left_on='target_id', right_on='transcript_id', how='inner')
    
    # Explode KEGG_ko column
    # First convert the string to a list
    data['annotation'] = data['annotation'].str.split(',')
    data = data.explode('annotation')
    
    # Run MWU
    res = enrich_mwu(data[['transcript_id', 'correlation']], data[['transcript_id', 'annotation']], min_size=5)

    # Add to results
    res['Cluster'] = cluster
    results = pd.concat([results, res])

# Save significant results
results = results[(results['q_val'] < 0.10) & (results['estimate'] > 0.5)]
results.to_csv('data/analysis/phaeocystis_O2_clusters_KEGG_MWU_results.csv', index=False)