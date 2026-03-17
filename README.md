# Transcriptome-resolved diel metabolism of a coastal *Phaeocystis* bloom links cellular physiology to oxygen production

## Description

This repository contains the code associated with the scientific manuscript titled 'Transcriptome-resolved diel metabolism of a coastal *Phaeocystis* bloom links cellular physiology to oxygen production'. The study leverages a multi-disciplinary approach to unravel the complexities of a late-stage *Phaeocystis globosa* bloom in the Southern North Sea. The study integrates biogeochemical data, metatranscriptomics, carbohydrate quantification, photophysiology, flowcam, and zooscan data to provide a comprehensive understanding of the bloom dynamics. The code provided in this repository supplements the manuscript by providing all the necessary steps to replicate the study's findings and analyses.

## Table of Contents

1. [Introduction](#introduction)
2. [Project structure](#Project-structure)
3. [Bioinformatics](#Bioinformatics)
4. [Analysis](#analyses)
5. [License](#license)

## Introduction
To the oceanographer, bloom situations present an opportunity to quantify the effect the blooming micro-algae have on the local biogeochemistry of the marine environment. *Phaeocystis globosa* is a cosmopolitan prymnesiophyte notorious for producing excessive amounts of marine gel during blooms. Short-lived *P. globosa* blooms heavily impact their environment, not only due to the formation of foam upon bloom termination or occasional toxin production, but also by high rates of primary production and nutrient drawdown. Through the combination of continuous biogeochemical measurements and hourly metatranscriptomics we assess the link between gene expression and in situ oxygen production and dissolved inorganic carbon (DIC) drawdown during a *P. globosa* bloom across a diel timeframe.

## Project structure

```bash
├── config.yaml
├── data
│   ├── analysis
│       ├── ...
│   └── raw
│       ├── ...
│
├── hpc_config
│   ├── cluster.yaml
│   └── config.yaml
├── README.md
├── figures
├── samples.csv
├── rules
│   ├── quality_control.smk
│   ├── preprocessing.smk
│   ├── assembly.smk
│   ├── cluster_assemblies.smk
│   ├── protein_prediction.smk
│   ├── quantification.smk
│   ├── spike_quantification.smk
│   ├── taxonomic_annotation.smk
│   ├── functional_annotation.smk
│   ├── phaglo1_mapping.smk
│   └── phaglo1_annotation.smk
├── scripts
│   ├── add_resequencing_files.sh
│   ├── biogeochemistry.R
│   ├── build_Phaglo1_gene_function_table.py
│   ├── build_tx2gene_from_gff.py
│   ├── combine_runs.sh
│   ├── environmental_analysis.ipynb
│   ├── ERCC_normalisation.ipynb
│   ├── flowcam.ipynb
│   ├── labstaf_processing.R
│   ├── labSTAF.R
│   ├── map.R
│   ├── mtx_taxonomy.ipynb
│   ├── phaeocystis_pathway_analysis.ipynb
│   ├── phaglo1_analysis.ipynb
│   ├── phaglo1_analysis.R
│   ├── primary_production_correlation.ipynb
│   ├── run_kallisto_merge.py
│   ├── satellite_chl_a.ipynb
│   ├── submit_merge_kallisto_pbs.sh
│   ├── submit_snakemake_pbs.sh
│   ├── sum_Phaglo1_transcripts_to_genes.py
│   ├── TEP_analysis.ipynb
│   ├── utils.py
│   ├── wgcna_enrichment_gene_functions.ipynb
│   └── zooscan.ipynb
└── Snakefile
```

## Bioinformatics
First, data from the two repeated sequencing runs is combined dynamically with this [script](scripts/combine_runs.sh). Then, Snakemake is deployed to process the metatranscriptomic data using this submission script that launches the [Snakefile](Snakefile). Snakemake is configured to run the following steps:
1. Quality control of the raw reads using FastQC and MultiQC.
2. Trimming of the raw reads using Trimmomatic.
3. rRNA removal using RiboDetector.
4. Assembly of the reads per sample using rnaSPAdes.
5. Clustering of the assembled transcripts using MMseqs2, generating the de novo metatranscriptome.
6. Prediction of the open reading frames using TransDecoder, and using the longest ones to translate the assembled transcripts into proteins.
7. Annotation of the metatranscripotome using the EUKProt reference database.
8. Functional annotation of the metatranscriptome using the eggNOG database and eggnog-mapper.
9. Quantification of the metatranscriptome using Kallisto.
10. Quantification of the ERCC Spike-ins using Kallisto and BBMap.
11. Mapping of the metatranscriptomic reads to the [Phaglo1 reference genome](https://phycocosm.jgi.doe.gov/Phaglo1/) using Salmon.
12. Gathering the Phaglo1 annotations and extending them with dbCAN annotations.

The resulting kallisto quantification files are merged using this [script](scripts/run_kallisto_merge.py). This generates a count.csv and tpm.csv file that are used in downstream analyses.

## Analyses
First, the environmental data is analysed in [this notebook](scripts/environmental_analysis.ipynb). In this notebook we integrate data from the nutrient analysis, tidal dynamics, data from the CTD casts, pull additional data using the [BPNSdata package](https://github.com/lifewatch/bpnsdata). Then we generate depth profiles of the CTD casts and T/S diagrams. This notebook generates the samples_env.csv file which is used in downstream analyses. A map of the sampling regions can be generated using the [map](scripts/map.R) script.

Underway data is processed in the [biogeochemistry script](scripts/biogeochemistry.r). This analysis calculates the oxygen saturation (O2') and DIC, and models and visualizes the diel patterns in these parameters.

The satellite-derived chlorophyll a concentrations are obtained and visualized [here](scripts/satellite_chl_a.ipynb).

[Sequencing QC and TPL calculation](scripts/ERCC_normalisation.ipynb) is done before processing the metatranscriptomic data. Relative and absolute taxonomic abundance plots are generated from the metatranscriptomic data in [this notebook](scripts/mtx_taxonomy.ipynb).

Flowcam data analysis is done [here](scripts/flowcam.ipynb). ZooScan data analysis is done [in this notebook](scripts/zooscan.ipynb).

The general patterns in the reads mapped to the Phaglo1 reference are analysed [here](scripts/phaglo1_analysis.ipynb). This script also generates a gene x TP table. The WGCNA on cyclically varying genes is performed in this [script](scripts/phaglo1_analysis.R). This is followed by [functional enrichment analyses](scripts/wgcna_enrichment_gene_functions.ipynb) on the obtained WGCNA clusters. The overall metabolic activity of *Phaeocystis globosa*, such as expression per KOG functional category, in station 130 and 51 is visualized [here](scripts/phaeocystis_pathway_analysis.ipynb).

[Photophysiology](scripts/labSTAF.R) analysis, using manually extracted LabSTAF data (all_Station_NF_final.csv) or [automatically extracted values](scripts/labstaf_processing.R). Correlations between environmental parameters, photophysiology, and gene expression are explored [here](scripts/primary_production_correlation.ipynb).

TEP analysis is done [here](scripts/TEP_analysis.ipynb).

## License
This code is licensed under the **Creative Commons Attribution 4.0 International (CC-BY 4.0)** license. See the [LICENSE](LICENSE) file for details.

### Citation
If we've inspired your analysis with this project, give us a shout out! You can cite us as follows:

Perneel & Dujardin, et al. "Transcriptome-resolved diel metabolism of a coastal *Phaeocystis* bloom links cellular physiology to oxygen production". [Journal Name], [2026]. DOI: [DOI]