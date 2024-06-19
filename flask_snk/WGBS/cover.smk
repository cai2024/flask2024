


rule qualimap_bamqc:
    input: 
        de_bam="{output_dir}/de_bam/{dup_soft}/sort_{{sample}}_dedup.bam".format(output_dir=config['output_dir'], dup_soft=config['dup_soft']),
    output: 
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/{{sample}}".format(output_dir=config['output_dir']),
    threads: 10
    log:
        "{output_dir}/logs/qualimap/{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.de_bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size=10G ) > {log} 2>&1
        """


