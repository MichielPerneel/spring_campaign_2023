import os
OUT = config['output_dir']
LOG = config['log_dir']

PHAGLO_GFF   = config['phaglo_gff']
PROT_FASTA   = config['phaglo_proteins']
PROT_TX_MAP  = config['phaglo_protein_transcript_map']
KEGG_TAB     = config['phaglo_kegg_tab']
GO_TAB       = config['phaglo_go_tab']
IPR_TAB      = config['phaglo_ipr_tab']
KOG_TAB      = config['phaglo_kog_tab']
SIGP_TAB     = config['phaglo_signalp_tab']
ENABLE_DBCAN = bool(config.get('enable_dbcan', False))

# (No tx2gene rule here; we depend on the one built in mapping.smk)
TX2GENE = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'tx2gene.tsv')

rule phaglo1_dbcan:
    """
    Optional CAZy annotation using dbCAN2; enable with config.enable_dbcan: true
    """
    input:
        proteins = PROT_FASTA
    output:
        overview = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'dbcan', 'overview.txt')
    params:
        outdir = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'dbcan'),
        dbdir  = config.get('dbcan_db_dir', os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'dbcan', 'db')),
        enable = "true" if ENABLE_DBCAN else "false"
    threads: 8
    conda: os.path.join(config['conda_environments'], 'dbcan.yaml')
    shell: r"""
    set -euo pipefail
    mkdir -p "{params.outdir}"

    # If disabled, write a stub and exit
    if [ "{params.enable}" = "false" ]; then
      printf "protein_id\tCAZy\n" > "{output.overview}"
      exit 0
    fi

    # Decompress proteins for dbCAN
    zcat '{input.proteins}' > "{params.outdir}/phaglo1_prots.fa"

    # Run dbCAN (v5.1.2)
    run_dbcan CAZyme_annotation \
      --mode protein \
      --input_raw_data "{params.outdir}/phaglo1_prots.fa" \
      --output_dir "{params.outdir}" \
      --db_dir "{params.dbdir}" \
      --methods hmm --methods diamond --methods dbCANsub \
      --threads {threads} \
      1> "{params.outdir}/run_dbcan.stdout" 2> "{params.outdir}/run_dbcan.stderr"

    # Keep downstream happy (your pipeline expects overview.txt)
    if [ -s "{params.outdir}/overview.tsv" ]; then
      cp "{params.outdir}/overview.tsv" "{output.overview}"
    else
      printf "protein_id\tCAZy\n" > "{output.overview}"
    fi
    """

rule phaglo1_build_gene_functions:
    """
    Join KEGG/GO/InterPro/KOG/SignalP (+optional dbCAN) at protein level,
    map to gene_id via protein->transcript->gene, and aggregate per gene.
    """
    input:
        tx2gene = TX2GENE,
        prot_tx = PROT_TX_MAP,
        kegg    = KEGG_TAB,
        go      = GO_TAB,
        ipr     = IPR_TAB,
        kog     = KOG_TAB,
        sigp    = SIGP_TAB,
        dbcan   = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'dbcan', 'overview.txt')
    output:
        func = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'gene_functions.tsv')
    log:
        os.path.join(LOG, 'phaglo1_mapping', 'build_gene_functions.log')
    conda: os.path.join(config['conda_environments'], 'salmon.yaml')
    threads: 1
    script: os.path.join(config['home_dir'], "scripts/build_Phaglo1_gene_function_table.py")