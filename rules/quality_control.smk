rule fastqc:
        input:
                os.path.join(config['raw_reads'], '{sample}_{num}.fastq.gz')
        output:
                html = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_{num}_fastqc.html'),
                zip = os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_{num}_fastqc.zip')
        log:
                os.path.join(config['output_dir'],'logs', 'fastqc', '{sample}_{num}.log')
        wrapper:
                '0.27.1/bio/fastqc'

rule multiqc:
    input:
        qc_files = expand(os.path.join(config['output_dir'], 'quality_control', 'fastqc', '{sample}_{num}_fastqc.{ext}'), sample=samples, num=[1, 2], ext=['html', 'zip'])
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

