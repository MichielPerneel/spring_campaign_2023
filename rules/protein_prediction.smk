rule protein_prediction:
        input: os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'metatranscriptome.fasta')
        output:
                bed=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.bed'),
                cds=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.cds'),
                gff3=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.gff3'),
                pep=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}', 'metatranscriptome.pep')
        params:
                longORFs_dir=os.path.join(config['output_dir'], 'assembly', 'protein', '{station}'),
                bed="metatranscriptome.fasta.transdecoder.bed",
                cds="metatranscriptome.fasta.transdecoder.cds",
                gff3="metatranscriptome.fasta.transdecoder.gff3",
                pep="metatranscriptome.fasta.transdecoder.pep",
                env = lambda wildcards: "module load TransDecoder" if config['environment_system'] == "modules" else ""
        log: os.path.join(config['log_dir'], 'protein_prediction_{station}.log')
        conda:
                os.path.join(config['conda_environments'], "transdecoder.yaml") if config['environment_system'] == "conda" else None
        shell: '''
                {params.env}
                cd {params.longORFs_dir}
                
                TransDecoder.LongOrfs -t {input} \
                        --output_dir {params.longORFs_dir}
                
                TransDecoder.Predict -t {input} \
                        --output_dir {params.longORFs_dir} --single_best_only
                
                mv {params.longORFs_dir}/{params.bed} {output.bed}
                mv {params.longORFs_dir}/{params.cds} {output.cds}
                mv {params.longORFs_dir}/{params.gff3} {output.gff3}
                mv {params.longORFs_dir}/{params.pep} {output.pep}
        '''