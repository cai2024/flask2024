rule picard_duplicate:
    input:
        bam="{output_dir}/bam/{align_soft}/sort_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"], align_soft=config['align_soft']),
    output:
        de_bam="{output_dir}/de_bam/picard/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
        sort_bam="{output_dir}/de_bam/picard/sort_{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/duplication/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    conda:
        "flask2024",
    params:
        temp_bam="{output_dir}/de_bam/picard/temp_{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
        tmp_report="{output_dir}/de_bam/{{sample}}.txt".format(output_dir=config['output_dir']),
        other=config['dup_params']
    shell:
        """
        ({java} -jar {picard} MarkDuplicates --REMOVE_DUPLICATES {params.other} -I {input.bam} -O {params.temp_bam} -M {params.tmp_report}) 2> {log}
        samtools sort -@ {threads}  {params.temp_bam}  > {output.sort_bam}    
        samtools index {output.sort_bam}
        samtools sort -n -@ {threads} {params.temp_bam}  > {output.de_bam}
        rm {params.temp_bam}
        """

