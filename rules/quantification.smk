rule kallisto_index:
    input:
        fasta = os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.fasta'),
    output:
        index = os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', '{station}', 'final_metatranscriptome.idx'),
    threads: 1
    wrapper: 'v1.15.0/bio/kallisto/index'

rule kallisto_quant:
    input:
        fastq = lambda wildcards: [
            os.path.join(config['scratch_dir'], 'cleanup', f'{wildcards.sample.strip()}_1.clean.fastq.gz').replace(" ", ""),
            os.path.join(config['scratch_dir'], 'cleanup', f'{wildcards.sample.strip()}_2.clean.fastq.gz').replace(" ", "")
        ],
        index = lambda wildcards: os.path.join(config['output_dir'], 'assembly', 'rnaSPAdes', f'{wildcards.sample.strip().split("_")[0]}', 'final_metatranscriptome.idx').replace(" ", ""),
    output: directory(os.path.join(config['output_dir'], 'quantification', '{sample}')),
    log: os.path.join(config['log_dir'], 'quantification', '{sample}.log'),
    threads: 1
    wrapper: 'v1.15.0/bio/kallisto/quant'