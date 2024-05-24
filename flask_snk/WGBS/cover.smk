rule qualimap_bamqc:
    input: 
        dedup_bam="{output_dir}/debam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'], dup_soft=config['dup_soft']),
    output: 
        "{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/{{sample}}".format(output_dir=config['output_dir']),
        java_mem_size="20G"
    threads: 10
    log:
        "{output_dir}/logs/qualimap/{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.dedup_bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size={params.java_mem_size}) > {log} 2>&1
        """


