rule kallisto_index_ERCC:
    input:
       fasta =  config['spikefile'],
    output:
        index = os.path.join(config['ERCC_folder'], 'ERCC92.idx'),
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/index"

rule kallisto_quant_ERCC:
    input:
        fastq = [
            os.path.join(config['scratch_dir'], 'cleanup', 
                        '{sample}_1.rRNA_removed.fastq.gz'),
            os.path.join(config['scratch_dir'], 'cleanup',
                        '{sample}_2.rRNA_removed.fastq.gz')],
        index = os.path.join(config['ERCC_folder'], 'ERCC92.idx'),
    output:
        directory(os.path.join(config['output_dir'],  'ERCC92', 'kallisto', '{sample}')),
    params:
        extra='',
    log:
        os.path.join(config['log_dir'], 'spike_quantification', '{sample}.log'),
    threads: 1
    wrapper:
        "v1.15.0/bio/kallisto/quant"