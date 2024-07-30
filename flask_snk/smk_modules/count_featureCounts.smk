if config["mode"] == "single":
    featureCounts_mode = ""
else:
    featureCounts_mode = "-p"



rule generate_expression_matrix:
    input:
        bam="{output_dir}/{{sample}}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"], align_soft=config['align_soft']),
    output:
        expression_matrix="{output_dir}/{{sample}}/count/featureCounts/{{sample}}_expression_matrix.txt".format(output_dir=config['output_dir']),
        log="{output_dir}/{{sample}}/logs/count/{{sample}}.log".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/{{sample}}/logs/count/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    conda:
        "flask2024",
    params:
        out_dir="{output_dir}/{{sample}}/count/featureCounts/".format(output_dir=config['output_dir']),
        other=config['count_params']
    shell:
        """
        (mkdir -p {params.out_dir} && \
        {featureCounts} {featureCounts_mode} -T {threads}  {params.other}  -a {gtf_path} -o {output.expression_matrix} {input.bam} ) >  {log} 2>&1 
        """



