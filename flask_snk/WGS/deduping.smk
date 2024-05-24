rule mark_duplicate:
    input:
        bam="{output_dir}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"], align_soft=config['align_soft']),
    output:
        dedup_bam="{output_dir}/debam/{align_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'], align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/duplication/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        tmp_report="{output_dir}/debam/{{sample}}.txt".format(output_dir=config['output_dir']),
        other=config['picard_params']
    shell:
        """
        (/usr/bin/java -jar {picard} MarkDuplicates --REMOVE_DUPLICATES {params.other} -I {input.bam} -O {output.dedup_bam} -M {params.tmp_report}) 2> {log}
        samtools index {output.dedup_bam}        
        """

