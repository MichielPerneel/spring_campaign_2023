import os

OUT = config['output_dir']
SCR = config['scratch_dir']
LOG = config['log_dir']
INCL_ORG = bool(config.get('include_organelles', False))

rule phaglo1_targets:
    """
    Build transcript target FASTA for Salmon: Phaglo1 transcripts (+ organelles if enabled).
    Currently we don't allow inclusion of organellar references, as these are genomic in the JGI reference.
    Assumes all inputs are .gz files.
    """
    input:
        phaglo_tx = config['phaglo_transcripts'],
        chl = config.get('phaglo_chloroplast', ''),
        mit = config.get('phaglo_mito', '')
    output:
        fa = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'phaglo1_targets.fa')
    params:
        incl_org = False   # keep organelles out of targets; they go to decoys
    threads: 1
    shell: r"""
        set -euo pipefail
        mkdir -p "{OUT}/phaglo1_mapping/ref"
        : > "{output.fa}"

        # nuclear transcripts
        zcat '{input.phaglo_tx}' >> "{output.fa}"

        # organelles NOT added to targets when params.incl_org is False
        if {params.incl_org}; then
          [ -n '{input.chl}' ] && zcat '{input.chl}' >> "{output.fa}"
          [ -n '{input.mit}' ] && zcat '{input.mit}' >> "{output.fa}"
        fi
    """

rule phaglo1_decoys:
    """
    Prepare decoy FASTA (repeat-masked nuclear genome + optional chl/mito genomes)
    and decoys.txt with unique names for Salmon selective alignment. 
    Using decoys avoids inflating transcript expression by reads that are actually genomic contamination
    (unspliced RNA, DNA contamination, repeats, organellar DNA reads)
    Chloroplast headers are prefixed with 'chl|' and mitochondrion with 'mit|' to avoid
    collisions with nuclear scaffold names.
    """
    input:
        masked = config['phaglo_masked_genome'],
        chl    = config.get('phaglo_chloroplast', ''),
        mit    = config.get('phaglo_mito', '')
    params:
        outdir = os.path.join(OUT, 'phaglo1_mapping', 'ref'),
    output:
        decoy_fa = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'phaglo1_decoys.fa'),
        decoys   = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'decoys.txt')
    threads: 1
    input:
        masked = config['phaglo_masked_genome'],
        chl    = config.get('phaglo_chloroplast', ''),
        mit    = config.get('phaglo_mito', '')
    params:
        outdir = os.path.join(OUT, 'phaglo1_mapping', 'ref'),
    output:
        decoy_fa = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'phaglo1_decoys.fa'),
        decoys   = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'decoys.txt')
    threads: 1
    shell: r"""
    set -euo pipefail
    mkdir -p "{params.outdir}"
    : > "{output.decoy_fa}"

    # 1) nuclear
    if [[ "{input.masked}" == *.gz ]]; then
      zcat "{input.masked}" >> "{output.decoy_fa}"
    else
      cat  "{input.masked}" >> "{output.decoy_fa}"
    fi

    # 2) chloroplast, prefix headers with chl|
    if [ -n "{input.chl}" ]; then
      if [[ "{input.chl}" == *.gz ]]; then
        zcat "{input.chl}" | sed 's/^>/>chl|/' >> "{output.decoy_fa}"
      else
        sed 's/^>/>chl|/' "{input.chl}" >> "{output.decoy_fa}"
      fi
    fi

    # 3) mitochondrion, prefix headers with mit|
    if [ -n "{input.mit}" ]; then
      if [[ "{input.mit}" == *.gz ]]; then
        zcat "{input.mit}" | sed 's/^>/>mit|/' >> "{output.decoy_fa}"
      else
        sed 's/^>/>mit|/' "{input.mit}" >> "{output.decoy_fa}"
      fi
    fi

    # Build decoys.txt
    sed -n 's/^>//p' "{output.decoy_fa}" > "{output.decoys}"

    # Sanity check for duplicate names
    dups=$(sed -n 's/^>//p' "{output.decoy_fa}" | sort | uniq -d | head -n 1 || true)
    if [ -n "$dups" ]; then
      echo "ERROR: duplicate decoy name detected after prefixing: $dups" >&2
      exit 1
    fi
    """

rule salmon_index_phaglo1:
    """
    Build Salmon selective-alignment index with decoys.
    """
    input:
        targets = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'phaglo1_targets.fa'),
        decoy_fa = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'phaglo1_decoys.fa'),
        decoys   = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'decoys.txt')
    output:
        idx = directory(os.path.join(OUT, 'phaglo1_mapping', 'salmon_index', 'phaglo1_sa')),
        combined = os.path.join(OUT, 'phaglo1_mapping', 'ref', 'targets_plus_decoys.fa')
    threads: 8
    conda: os.path.join(config['conda_environments'], 'salmon.yaml')
    shell: r"""
    set -euo pipefail
    cat "{input.targets}" "{input.decoy_fa}" > "{output.combined}"
    salmon index -t "{output.combined}" -d "{input.decoys}" -i "{output.idx}" -k 23
    """

rule salmon_quant_phaglo1:
    """
    Quantify each sample against Phaglo1 with decoy-aware Salmon.
    """
    input:
        idx = os.path.join(OUT, 'phaglo1_mapping', 'salmon_index', 'phaglo1_sa'),
        r1  = lambda wc: os.path.join(SCR, 'cleanup', f"{wc.sample.strip()}_1.clean.fastq.gz").replace(" ", ""),
        r2  = lambda wc: os.path.join(SCR, 'cleanup', f"{wc.sample.strip()}_2.clean.fastq.gz").replace(" ", "")
    output:
        sf = os.path.join(OUT, 'phaglo1_mapping', 'salmon_quant', '{sample}', 'quant.sf')
    params:
        outdir = os.path.join(OUT, 'phaglo1_mapping', 'salmon_quant', '{sample}'),
        libtype = config.get('salmon_libtype', 'ISR')
    log:
        os.path.join(LOG, 'phaglo1_mapping', 'salmon', '{sample}.log')
    threads: 8
    conda: os.path.join(config['conda_environments'], 'salmon.yaml')
    shell: r"""
    set -euo pipefail
    mkdir -p "{params.outdir}" "$(dirname "{log}")"
    salmon quant -i "{input.idx}" \
      -l {params.libtype} \
      -1 "{input.r1}" -2 "{input.r2}" \
      --validateMappings --seqBias --gcBias --posBias \
      -p {threads} -o "{params.outdir}" &> "{log}"
    """

rule phaglo1_tx2gene:
    """
    Build tx2gene table from the JGI GFF (transcript/mRNA -> gene).
    """
    input:
        gff = config['phaglo_gff']
    output:
        tx2gene = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'tx2gene.tsv')
    conda: os.path.join(config['conda_environments'], 'salmon.yaml')
    log:    os.path.join(LOG, 'phaglo1_mapping', 'tx2gene.log')
    threads: 1
    script: os.path.join(config['home_dir'], "scripts/build_tx2gene_from_gff.py")

rule phaglo1_merge_gene_level:
    """
    Sum Salmon transcript TPM/NumReads to gene-level; merge across all samples.
    """
    input:
        tx2gene = os.path.join(OUT, 'phaglo1_mapping', 'annotation', 'tx2gene.tsv'),
        sfs = expand(os.path.join(OUT, 'phaglo1_mapping', 'salmon_quant', '{sample}', 'quant.sf'), sample=samples)
    output:
        tpm = os.path.join(OUT, 'phaglo1_mapping', 'gene_expression', 'gene_tpm.tsv'),
        numreads = os.path.join(OUT, 'phaglo1_mapping', 'gene_expression', 'gene_numreads.tsv')
    conda: os.path.join(config['conda_environments'], 'salmon.yaml')
    threads: 1
    script: os.path.join(config['home_dir'], "scripts/sum_Phaglo1_transcripts_to_genes.py")