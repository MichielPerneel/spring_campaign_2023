# Metabolic Insights into Diel Carbon Cycling and Oxygen Production during a *Phaeocystis globosa* bloom

## Description

This repository contains the code associated with the scientific manuscript titled 'Metabolic Insights into Diel Carbon Cycling and Oxygen Production during a *Phaeocystis globosa* bloom'. The study leverages a multi-disciplinary approach to unravel the complexities of a late-stage *Phaeocystis globosa* bloom in the North Sea. The study integrates biogeochemical data, metatranscriptomics, carbohydrate quantification, photophysiology, flowcam, and zooscan data to provide a comprehensive understanding of the bloom dynamics. The code provided in this repository is aimed at supplementing the manuscript by providing the necessary material to replicate the study's findings and analyses.

## Table of Contents

1. [Introduction](#introduction)
2. [Project structure](#Project-structure)
3. [Bioinformatics](#Bioinformatics)
4. [Analysis](#analyses)
5. [License](#license)

## Introduction
To the oceanographer, bloom situations present an opportunity to quantify the effect the blooming micro-algae have on the local biogeochemistry of the marine environment. Phaeocystis globosa is a cosmopolitan prymnesiophyte notorious for producing excessive amounts of marine gel during blooms. Short-lived P. globosa blooms heavily impact their environment, not only due to the formation of foam upon bloom termination or occasional toxin production, but also by high rates of primary production and nutrient drawdown. Through the combination of biogeochemical measurements and metatranscriptomics we are able to observe in situ O2 production and dissolved inorganic carbon (DIC) drawdown by P. globosa across a diel timeframe. Using synchronized hourly sampling, we are able to propose a set of metabolic modules and potential genetic markers whose expression closely follows biogeochemical changes. The results from this case study contribute towards describing physiological and biogeochemical rates from omics data.

## Project structure

```bash
├── config.yaml
├── data
│   ├── analysis
│       ├── 
│   └── raw
│       ├── 
│
├── hpc_config
│   ├── cluster.yaml
│   └── config.yaml
├── README.md
├── figures
├── samples.csv
├── scripts
│   ├── TEP_analysis.ipynb
│   ├── carbonate_chemistry.R
│   ├── combine_runs.sh
│   ├── environmental_analysis.ipynb
│   ├── ERCC_normalisation.ipynb
│   ├── flowcam.ipynb
│   ├── get_marker_transcripts.py
│   ├── labstaf_processing.R
│   ├── map.R
│   ├── mtx_taxonomy.ipynb
│   ├── oxygen_light_tides.R
│   ├── phaeocystis_cluster_enrichment_visualization.ipynb
│   ├── photophysiology.ipynb
│   ├── run_kallisto_merge.py
│   ├── submit_cluster_MWU.pbs
│   ├── submit_marker_detection.pbs
│   ├── submit_phaeo_clustering.pbs
│   ├── submit_snakemake_pbs.sh
│   ├── transcript_markers.ipynb
│   ├── zooscan.ipynb
│   └── mbcluster_phaeocystis.R
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

The resulting kallisto quantification files are merged using this [script](scripts/run_kallisto_merge.py). This generates a count.csv and tpm.csv file that are used in downstream analyses.

## Analyses
First, the environmental data is analysed in [this notebook](scripts/environmental_analysis.ipynb). In this notebook we integrate data from the nutrient analysis, 
tidal dynamics, data from the CTD casts, pull additional data using the [BPNSdata package](https://github.com/lifewatch/bpnsdata). Then we generate depth profiles of the CTD casts and T/S diagrams. This notebook generates the samples_env.csv file which is used in downstream analyses. A map of the sampling regions can be generated using the [map](scripts/map.R) script. 

Underway data is processed in the[ Carbonate chemistry script](scripts/carbonate_chemistry.R). This analysis extracts values from the fitted smoother functions as well, which are used downstream to relate the omics findings with the biogeochemical patterns.

Now, we can run an [analysis](scripts/oxygen_light_tides.R) of the light and tides, producing a figure that shows how dissolved oxygen concentration (from the smoothers generated before) changes with light and tide.

[Sequencing QC and TPL calculation](scripts/ERCC_normalisation.ipynb) is done before processing the metatranscriptomic data. Relative and absolute taxonomic abundance plots are generated from the metatranscriptomic data in [this notebook](scripts/mtx_taxonomy.ipynb).

Flowcam data analysis is done [here](scripts/flowcam.ipynb). ZooScan data analysis is done [in this notebook](scripts/zooscan.ipynb).

[Photophysiology](scripts/labSTAF.R) analysis, using manually extracted LabSTAF data or [automatically extracted values](scripts/labstaf_processing.R).

The metabolic activity of *Phaeocystis globosa* in station 130 and 51 is visualized [here](scripts/phaeocystis_pathway_analysis.ipynb).

Marker gene detection and clustering of expressed genes in high/low dissolved oxygen conditions is done in two approaches. Identifying sign. potential marker genes is done by submitting the [get_marker_transcripts.py](scripts/get_marker_transcripts.py) script to the HPC using [this submission script](scripts/submit_marker_detection.pbs). Clustering genes according to high/low dissolved oxygen conditions is done by submitting the [mbcluster_phaeocystis.R](scripts/mbcluster_phaeocystis.R) script to the HPC using [this submission script](scripts/submit_phaeo_clustering.pbs). In [this notebook](scripts/marker_gene_visualization.R) we visualize marker genes for high and low oxygen conditions. Information is extracted from the results in this [notebook](scripts/transcript_markers.ipynb). During that analysis, a correlation of the marker genes with the dissolved oxygen concentration is done, which is necessary for doing the MWU ranked tests to identify metabolic pathways that are differentially expressed in high and low oxygen conditions. If you have that file, you can run the [MWU ranked test](scripts/cluster_enrichment_MWU) to identify the metabolic pathways that are differentially expressed in high and low oxygen conditions, by submitting this [script](scripts/submit_cluster_MWU.pbs) on the HPC. MWU ranked test results are visualized in [this notebook](scripts/phaeocystis_cluster_enrichment_visualization.ipynb).

TEP analysis is done [here](scripts/TEP_analysis.ipynb).


## License
This code is licensed under the **Creative Commons Attribution 4.0 International (CC-BY 4.0)** license. See the [LICENSE](LICENSE) file for details.

### Citation
If we've inspired your analysis with our project, give us a shout out! You can cite us as follows:

**Perneel & Dujardin, et al. "Metabolic Insights into Diel Carbon Cycling and Oxygen Production during a *Phaeocystis globosa* bloom." [Journal Name], [2025]. DOI: [DOI]**