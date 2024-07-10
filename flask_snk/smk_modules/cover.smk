rule qualimap_bamqc:
    input: 
        bam="{output_dir}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'], align_soft=config['align_soft'], mode=config["mode"]),
    output: 
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/{{sample}}".format(output_dir=config['output_dir']),
    threads: 10
    log:
        "{output_dir}/logs/qualimap/{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size=10G ) > {log} 2>&1
        """

