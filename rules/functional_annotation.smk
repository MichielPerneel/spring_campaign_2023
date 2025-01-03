rule functional_annotation:
    input:
        pep=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep'),
    output:
        annotated=os.path.join(config['output_dir'], 'annotation', 'functional_eggnog', '{station}', 'functional_annotation.emapper.annotations'),
    params:
        eggnog_ref_dir=config['eggnog_ref_dir'],
        eggnog_out_dir=os.path.join(config['output_dir'], 'annotation', 'functional_eggnog', '{station}'),
        env = lambda wildcards: "module load eggnog-mapper" if config['environment_system'] == "modules" else ""
    threads: 48
    resources:
        mem_mb=100000
    log: os.path.join(config['log_dir'], 'functional_annotation_{station}.log')
    conda:
        os.path.join(config['conda_environments'], "eggnog-mapper.yaml") if config['environment_system'] == "conda" else None
    shell: '''
    {params.env}

    unset OMP_PROC_BIND

    mkdir -p {params.eggnog_out_dir}

    emapper.py -m mmseqs --no_annot --no_file_comments --cpu {threads} -i {input.pep} \
        -o eggNOG --output_dir {params.eggnog_out_dir} \
        --data_dir {params.eggnog_ref_dir} --override

    emapper.py --annotate_hits_table {params.eggnog_out_dir}/eggNOG.emapper.seed_orthologs \
        --no_file_comments --data_dir {params.eggnog_ref_dir} -o functional_annotation \
        --cpu {threads} --dbmem --output_dir {params.eggnog_out_dir} --override

    mv {params.eggnog_out_dir}/functional_annotation.emapper.annotations {output.annotated}
    '''
