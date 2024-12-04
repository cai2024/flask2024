

rule clipOverlap:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bam = "{output_dir}/{sample}/{step}_clipOverlap/clipOverlap_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_clipOverlap/{sample}.log"
    threads: 4
    conda:
        "flask2024",
    params:
        bam = config["bamUtil"],
        tmp_bam = "{output_dir}/{sample}/{step}_clipOverlap/clipOverlap_{sample}_temp.bam",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        final_log="{output_dir}/logs/{sample}.log",
    shell:
        """
        samtools sort -n {input.bam} -@ {threads} -o {params.tmp_bam}
        ({params.bam} clipOverlap --readName --in {params.tmp_bam} --out - | \
        samtools sort -@ {threads}  -o {output.bam} -O bam - ) 2> {log}
        samtools index {output.bam}
        rm {params.tmp_bam} 
        """

