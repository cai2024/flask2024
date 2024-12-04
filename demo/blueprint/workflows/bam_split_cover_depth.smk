
rule bam_split_cover_depth:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"]
    output:
        txt = "{output_dir}/{sample}/{step}_bam_split_cover_depth/bam_split_cover_depth_{sample}.txt"
    params:
        output="{output_dir}/{sample}/{step}_bam_split_cover_depth/bam_split_cover_depth_{sample}",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        qualimap = config["qualimap"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bam_split_cover_depth/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/plot_bam_split_cover_depth.py --bam {input.bam} --output {params.output} --qualimap {params.qualimap} {params.other} > {log} 2>&1

        """
