configfile: "config.yaml"

import os
import sys
import pandas as pd
# Add 'scripts' and its subdirectories to Python's module search path to allow importing custom modules. 
sys.path.insert(1, 'scripts')

from utils import load_config, extract_samplenames, check_config

# Load configuration.
config = load_config('config.yaml')

# Load samples.
samples = extract_samplenames(config["samples"])
# Load the sample information
samples_df = pd.read_csv("data/samples.csv", sep=';', decimal=',')
# Create a dictionary mapping from sample name to station
sample_to_station_map = pd.Series(samples_df.station.values,index=samples_df.sample_name).to_dict()
stations = set(sample_to_station_map.values())
stations = set(str(station) for station in stations)

# Check setup
check_config(config)

# Include rules from separate modules.
include: "rules/quality_control.smk"
include: "rules/preprocessing.smk"
include: "rules/assembly.smk"
include: "rules/cluster_assemblies.smk"
include: "rules/protein_prediction.smk"
include: "rules/quantification.smk"
include: "rules/spike_quantification.smk"
include: "rules/taxonomic_annotation.smk"
include: "rules/functional_annotation.smk"

ruleorder: bbmap_combined > combined_sample_rRNA_cleanup > trimmomatic_combined > bbmap > sample_rRNA_cleanup > trimmomatic

# Run rules.
rule all:
       input:
              #expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.html'),
              #       sample=samples, num=[1, 2]),
              #expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.zip'),
              #       sample=samples, num=[1, 2]),
              #os.path.join(config['output_dir'], 'quality_control', 'multiqc_data'),
              expand(os.path.join(config['output_dir'], "assembly", "rnaSPAdes", "{sample}", "transcripts.fasta"), sample=samples),
              expand(os.path.join(config['ERCC_folder'], 'ERCC92.idx')),
              expand(os.path.join(config['output_dir'], 'ERCC92', 'kallisto', '{sample}'), sample=samples),
              expand(os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.fasta'), station=stations),
              expand(os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep'), station=stations),
              expand(os.path.join(config['output_dir'], 'annotation', 'functional_eggnog', '{station}', 'functional_annotation.emapper.annotations'), station=stations),
              expand(os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot', '{station}', 'eukprot_annotation.m8'), station=stations),
              expand(os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.idx'), station=stations),
              expand(os.path.join(config['output_dir'], 'quantification', '{sample}'), sample=samples)