rule mark_duplicate:
    input:
        bam="{output_dir}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"], align_soft=config['align_soft']),
    output:
        de_bam="{output_dir}/de_bam/picard/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/duplication/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    conda:
        "flask2024",
    params:
        tmp_report="{output_dir}/de_bam/picard/{{sample}}.txt".format(output_dir=config['output_dir']),
        other=config['dup_params']
    shell:
        """
        ({java} -jar {picard} MarkDuplicates --REMOVE_DUPLICATES {params.other} -I {input.bam} -O {output.de_bam} -M {params.tmp_report}) > {log} 2>&1
        samtools index {output.de_bam}        
        """

