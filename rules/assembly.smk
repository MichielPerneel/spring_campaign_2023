rule rnaSPAdes:
        input:
                r1 = os.path.join(config['scratch_dir'], 'cleanup',
                         'combined', '{sample}_1.clean.fastq.gz'),
                r2 = os.path.join(config['scratch_dir'], 'cleanup',
                         'combined', '{sample}_2.clean.fastq.gz') 
        output:
                dir = directory(os.path.join(config['output_dir'], 
                                'assembly', 'rnaSPAdes', '{sample}')),
                transcripts = os.path.join(config['output_dir'], 'assembly',
                                'rnaSPAdes', '{sample}', 'transcripts.fasta'),
        log:
                err = os.path.join(config['output_dir'], 'logs', 'rnaspades', '{sample}_rnaspades.err'),
                out = os.path.join(config['output_dir'], 'logs', 'rnaspades', '{sample}_rnaspades.out')
        params:
                k = config['k_values']
        shell:'''
        module load SPAdes
        spades.py --rna \
                -1 {input.r1} \
                -2 {input.r2} \
                -o {output.dir} \
                -k {params.k} \
                2> {log.err} 1> {log.out}
        '''