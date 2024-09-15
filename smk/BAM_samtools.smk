rule mark_primary:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bam = "{output_dir}/{sample}/{step}_primary_filtering/primary_filtering_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_primary_filtering/{sample}.log"
    threads: 4
    conda:
        "flask2024",
    params:
        other = lambda wildcards: config[f"{wildcards.step}_params"],
    shell:
        """
        samtools view -bh -F 256 -F 2048 {params.other} -@ {threads} {input.bam} > {output.bam}
        """

