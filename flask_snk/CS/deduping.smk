rule mark_duplicate: 
    input:
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']) ,
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        temp_bam=temp("{output_dir}/debam/{align_soft}/de_temp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft'])),
        bam="{output_dir}/debam/{align_soft}/de_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft'])
    log:
        "{output_dir}/logs/dedup/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 10
    shell:
        """
        (python3 {py_ref}/remove_dup.py --input {input.bam} --output {output.temp_bam} --fq1 {input.fq1} --fq2 {input.fq2} ) 2> {log}
        samtools sort -@ {threads} {output.temp_bam} -o {output.bam}
        samtools index {output.bam}
        """

