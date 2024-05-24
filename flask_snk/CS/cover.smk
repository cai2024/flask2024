rule qualimap_bamqc:
    input: 
        dedup_bam="{output_dir}/debam/{align_soft}/de_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    output: 
        "{output_dir}/cover/de_{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/de_{{sample}}".format(output_dir=config['output_dir']),
        java_mem_size="20G"
    threads: 10
    log:
        "{output_dir}/logs/qualimap/de_{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.dedup_bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size={params.java_mem_size}) > {log} 2>&1
        """

