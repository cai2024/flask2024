rule picard_duplicate:
    input:
        bam="{output_dir}/bam/{align_soft}/sort_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"], align_soft=config['align_soft']),
    output:
        dedup_bam="{output_dir}/debam/picard/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/duplication/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        temp_bam="{output_dir}/debam/picard/temp_{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
        tmp_report="{output_dir}/debam/{{sample}}.txt".format(output_dir=config['output_dir']),
        other=config['dup_params']
    shell:
        """
        (/usr/bin/java -jar {picard} MarkDuplicates --REMOVE_DUPLICATES {params.other} -I {input.bam} -O {params.temp_bam} -M {params.tmp_report}) 2> {log}
        samtools sort {params.temp_bam}  > {output.dedup_bam}    
        samtools index {output.dedup_bam}
        rm {params.temp_bam}
        """

