rule cluster_assemblies:
        input: 
                lambda wildcards: expand(os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{sample}', 'transcripts.fasta'), sample=[sample for sample, station in sample_to_station_map.items() if str(station) == str(wildcards.station)])
        output:
                metatranscriptome=os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'metatranscriptome.fasta'),
        params:
                percent_id = "0.98",
                mmseqs_clu=os.path.join(config['scratch_dir'], 'mmseqs_clustering', '{station}', 'clusterDB'),
                tmp=os.path.join(config['scratch_dir'], 'mmseqs_clustering', '{station}', 'scratch'),
                env = lambda wildcards: "module load MMseqs2" if config['environment_system'] == "modules" else "",
        log: os.path.join(config['log_dir'], 'cluster_assemblies_{station}.log')
        conda: os.path.join(config['conda_environments'], "mmseqs2.yaml") if config['environment_system'] == "conda" else None
        shell: '''
        {params.env}

        unset OMP_PROC_BIND

        mmseqs easy-linclust {input} {params.mmseqs_clu} {params.tmp} --min-seq-id {params.percent_id} > {log} 2>&1

        mv {params.mmseqs_clu}_rep_seq.fasta {output.metatranscriptome}
        '''

rule rename_transcripts:
        input: os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'metatranscriptome.fasta')
        output: os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.fasta')
        params: env = lambda wildcards: "module load anvio" if config['environment_system'] == "modules" else "",
        log: os.path.join(config['log_dir'], 'rename_transcripts_{station}.log')
        shell: '''
        {params.env}
        anvi-script-reformat-fasta {input} -o {output} --simplify-names > {log} 2>&1
        '''