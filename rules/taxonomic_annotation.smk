rule create_eukprot_db:
    """
    This rule turns the EukProt sequences into a mmseqs DB
    """
    input: config['eukprot_ref_fa']
    output: os.path.join(config['eukprot_ref_dir'], 'eukprot')
    shell:'''
    unset OMP_PROC_BIND

    module load MMseqs2

    mmseqs createdb {input} {output}
    '''

rule create_metatranscriptome_db:
    """
    This rule turns the metatranscriptome sequences into a mmseqs DB for each station.
    """
    input: 
        os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep')
    output: 
        os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome_mmseqsDB')
    params:
        env = lambda wildcards: "module load MMseqs2" if config['environment_system'] == "modules" else ""
    conda:
        os.path.join(config['conda_environments'], "mmseqs2.yaml") if config['environment_system'] == "conda" else None
    shell: '''
    {params.env}
    unset OMP_PROC_BIND

    mmseqs createdb {input} {output}
    '''

rule eukprot_annotation:
    """
    This rule annotates the metatranscriptome sequences against the EukProt database for each station.
    """
    input:
        metatranscriptome_mmseqsDB=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome_mmseqsDB'),
        eukprot_mmseqsDB=os.path.join(config['eukprot_ref_dir'], 'eukprot')
    output:
        os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot', '{station}', 'eukprot_annotation.m8')
    params:
        eukprot_out_dir=os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukprot', '{station}'),
        tmp=config['scratch_dir'],
        env = lambda wildcards: "module load MMseqs2" if config['environment_system'] == "modules" else ""
    threads: 60
    resources:
        mem_mb=250000
    conda:
        os.path.join(config['conda_environments'], "mmseqs2.yaml") if config['environment_system'] == "conda" else None
    log: os.path.join(config['log_dir'], 'eukprot_annotation_{station}.log')
    shell: '''
    {params.env}
    unset OMP_PROC_BIND

    mkdir -p {params.eukprot_out_dir}

    # Query assembly against reference
    mmseqs search {input.metatranscriptome_mmseqsDB} {input.eukprot_mmseqsDB} {params.eukprot_out_dir}/resultDB {params.tmp} -s 6

    # Extract first hit
    mmseqs filterdb {params.eukprot_out_dir}/resultDB {params.eukprot_out_dir}/resultDB.firsthit --extract-lines 1
    mmseqs convertalis {input.metatranscriptome_mmseqsDB} {input.eukprot_mmseqsDB} {params.eukprot_out_dir}/resultDB.firsthit {output}
    '''

rule create_marferret_db:
    input:
        fasta=config['marferret_ref_fa'],
        taxonnodes=os.path.join(config['marferret_ref_dir'], 'ncbi/nodes.dmp'),
        taxonnames=os.path.join(config['marferret_ref_dir'], 'ncbi/names.dmp'),
        taxonmap=os.path.join(config['marferret_ref_dir'], 'MarFERReT.v1.taxonomies.tab.gz'),
    output:
        db=os.path.join(config['marferret_ref_dir'], 'dmnd/MarFERReT.v1.dmnd'),
    shell:'''
    module load DIAMOND/2.1.8-GCC-12.3.0
    diamond makedb --in {input.fasta} --db {output.db} --taxonnodes {input.taxonnodes} \
    --taxonnames {input.taxonnames} --taxonmap {input.taxonmap}
    '''

rule marferret_annotation:
    """
    This rule annotates the metatranscriptome sequences against the MarFERReT database for each station.
    """
    input:
        db=os.path.join(config['marferret_ref_dir'], 'dmnd/MarFERReT.v1.dmnd'),
        seqs=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep')
    output:
        annotated=os.path.join(config['output_dir'], 'annotation', 'taxonomy_marferret', '{station}', 'metatranscriptome_vs_MarFERReT_v1.lca.tab')
    params:
        evalue=config['evalue'],
        env = lambda wildcards: "module load DIAMOND/2.1.8-GCC-12.3.0" if config['environment_system'] == "modules" else ""
    threads: 60
    resources:
        mem_mb=200000
    conda:
        os.path.join(config['conda_environments'], "diamond.yaml") if config['environment_system'] == "conda" else None
    log: os.path.join(config['log_dir'], 'marferret_annotation_{station}.log')
    shell: '''
    {params.env}
    diamond blastp -b 150 -c 1 -d {input.db} -e {params.evalue} --top 10 -f 102 -q {input.seqs} -o {output.annotated} > {log} 2>&1
    '''

rule eukulele:
    input:
        pep = os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep')
    output:
        eukulele_done = os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukulele', '{station}', 'EUKulele_done.txt')
    params:
        sampledir = os.path.join(config['output_dir'], 'assembly', 'protein', '{station}'),
        eukulele_directory = os.path.join(config['output_dir'], 'annotation', 'taxonomy_eukulele', '{station}'),
        eukulele_reference_dir = config['eukulele_ref_dir'],
        env = lambda wildcards: "module load EUKulele" if config['environment_system'] == "modules" else ""
    conda: os.path.join(config['conda_environments'], "eukulele.yaml") if config['environment_system'] == "conda" else None
    shell: '''
    {params.env}
    EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir {params.eukulele_reference_dir} -o {params.eukulele_directory}
    touch {output.eukulele_done}
    '''