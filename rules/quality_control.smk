rule fastqc:
    input:
        os.path.join(config['raw_reads'], '{sample}_R{num}.fastq.gz')
    output:
        html = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.html'),
        zip = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.zip')
    log: 
        os.path.join(config['output_dir'], 'logs', 'fastqc', '{sample}_R{num}.log')
    shell: '''
    module load FastQC
    fastqc {input} -o {config[output_dir]}/quality_control/fastqc &> {log}
    '''

rule multiqc:
    input:
        qc_files = expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_R{num}_fastqc.{ext}'), sample=samples, num=[1, 2], ext=['html', 'zip'])
    output:
        report = os.path.join(config['output_dir'], 'quality_control', 'multiqc_report.html'),
        data_dir = directory(os.path.join(config['output_dir'], 'quality_control', 'multiqc_data'))
    log:
        os.path.join(config['output_dir'], 'logs', 'multiqc', 'multiqc.log')
    shell:
        """
        module load MultiQC
        multiqc {input} \
                --force \
                -o {output.data_dir} \
                -n {output.report} \
                &> {log}
        """