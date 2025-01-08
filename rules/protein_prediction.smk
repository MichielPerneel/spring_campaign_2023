rule protein_prediction:
        """
        This rule predicts proteins from the metatranscriptome as follows:
        1. Predicts long ORFs using TransDecoder.LongOrfs
        2. Searches for homologs using diamond blastp (optional)
        3. Predicts proteins using TransDecoder.Predict
        """
        input: os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.fasta')
        output:
                bed=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.bed'),
                cds=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.cds'),
                gff3=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.gff3'),
                pep=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep')
        params:
                longORFs_dir=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}'),
                bed="final_metatranscriptome.fasta.transdecoder.bed",
                cds="final_metatranscriptome.fasta.transdecoder.cds",
                gff3="final_metatranscriptome.fasta.transdecoder.gff3",
                pep="final_metatranscriptome.fasta.transdecoder.pep",
                env = lambda wildcards: """
                        module load TransDecoder
                        module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0
                        """ if config['environment_system'] == "modules" else ""
        log: os.path.join(config['log_dir'], 'protein_prediction_{station}.log')
        conda:
                os.path.join(config['conda_environments'], "transdecoder.yaml") if config['environment_system'] == "conda" else None
        shell:'''
        {params.env}

        cd {params.longORFs_dir}

        TransDecoder.LongOrfs -t {input} \
                --output_dir {params.longORFs_dir} > {log} 2>&1
        
        TransDecoder.Predict -t {input} \
                --output_dir {params.longORFs_dir} \
                --single_best_only >>{log} 2>&1

        mv {params.longORFs_dir}/{params.bed} {output.bed} >> {log} 2>&1
        mv {params.longORFs_dir}/{params.cds} {output.cds} >> {log} 2>&1
        mv {params.longORFs_dir}/{params.gff3} {output.gff3} >> {log} 2>&1
        mv {params.longORFs_dir}/{params.pep} {output.pep} >> {log} 2>&1
        '''