rule bowtie2_se:
    input:
        fq=lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam="{output_dir}/{sample}/{step}_bowtie2/bowtie2_{sample}.bam"
    params:
        index=lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bowtie2/" +
            config[f"{wildcards.step}_ref_release"]
        ),
        other=lambda wildcards: config[f"{wildcards.step}_params"],
        bowtie2 = config["bowtie2"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bowtie2/{sample}.log"
    threads: 10
    shell:
        """
        {params.bowtie2} {params.other} -x {params.index} -U {input.fq} -p {threads} | \
        samtools view -bS -F 256 -F 2048 - | samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}

        total_reads=$(zcat {input.fq} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(100 * $aligned_reads / $total_reads, '.2f'))" )      
        echo "Alignment Percentage $alignment_percentage%" >> {params.final_log}
        """
rule bowtie2_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bowtie2/bowtie2_{sample}.bam"
    params:
        index = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bowtie2/" +
            config[f"{wildcards.step}_ref_release"]
        ),
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bowtie2 = config["bowtie2"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bowtie2/{sample}.log"
    threads: 10
    shell:
        """
        {params.bowtie2} {params.other} -x {params.index} -1 {input.fq1} -2 {input.fq2} -p {threads} | \
        samtools view -bS -F 256 -F 2048 - | samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}

        total_reads=$(zcat {input.fq1} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(50 * $aligned_reads / $total_reads, '.2f'))" ) 
        echo "Alignment Percentage $alignment_percentage%" >> {params.final_log}
        """

rule sp_bowtie2_se:
    input:
        fq=lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam="{output_dir}/{sample}/{step}_bowtie2_spikein/bowtie2_spikein_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_bowtie2_spikein/{sample}.log"
    threads: 10
    conda:
        "flask2024",
    params:
        index = config["ref_path"] + "/spikein/bowtie2/spikein",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bowtie2 = config["bowtie2"],
    shell:
        """
        {params.bowtie2} {params.other} -x {params.index} -U {input.fq} -p {threads} | \
        samtools view -bS - | samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}
        """

rule sp_bowtie2_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bowtie2_spikein/bowtie2_spikein_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_bowtie2_spikein/{sample}.log"
    threads: 10
    conda:
        "flask2024"
    params:
        index = config["ref_path"] + "/spikein/bowtie2/spikein",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bowtie2 = config["bowtie2"]
    shell:
        """
        {params.bowtie2} {params.other} -x {params.index} -1 {input.fq1} -2 {input.fq2} -p {threads} | \
        samtools view -bS - | samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}
        """


