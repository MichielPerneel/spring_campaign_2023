# Unlocking the Secrets of Transparent Exopolymer Particles (TEP) Production in Coastal Planktonic Ecosystems: A Multi-disciplinary Approach

## Description

This repository contains the code associated with the scientific manuscript titled 'Unlocking the Secrets of Transparent Exopolymer Particles (TEP) Production in Coastal Planktonic Ecosystems: A Multi-disciplinary Approach'. The study leverages a multi-disciplinary approach to unravel the complexities of TEP production in coastal planktonic ecosystems.

## Table of Contents

1. [Introduction](#introduction)
2. [Project structure](#Project-structure)
4. [Contributing](#contributing)
5. [License](#license)

## Introduction

Detailed background and the objective of the study are provided in the manuscript. This code repository aims to supplement the manuscript by providing the necessary code to replicate the study's findings and analyses.

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
└── snakefile
```

## Contributing

1. Fork it (<https://github.com/yourusername/yourrepositoryname/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Create and activate the conda environment

Use the environment.yml file to create a dedicated conda environment. This file includes all the necessary dependencies for this project.

```bash
conda env create -f environment.yaml
```

Once the environment is created, you can activate it using:

```bash
conda activate spring_campaign_2023
```

If you need to updates to the environment, you can update the environment.yml file and run the following command:

```bash
conda env update -f environment.yaml  --prune
```

## Analyses
TEP analysis is done [here](scripts/TEP_analysis.ipynb).
Flowcam data analysis is done [here](scripts/flowcam.ipynb).

## License


## Acknowledgments
