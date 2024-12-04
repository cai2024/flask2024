rule samtools_merge:
    input:
        bam1 = lambda wildcards: config["steps"][wildcards.step]["input"]["bam1"],
        bam2 = lambda wildcards: config["steps"][wildcards.step]["input"]["bam2"],
    output:
        bam = "{output_dir}/{sample}/{step}_samtools_merge/samtools_merge_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_samtools_merge/{sample}.log"
    threads: 4
    conda:
        "flask2024",
    params:
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        final_log="{output_dir}/logs/{sample}.log",
    shell:
        """
        samtools  merge {output.bam} {input.bam1} {input.bam2} {params.other}
        samtools index {output.bam}    
   
        """

