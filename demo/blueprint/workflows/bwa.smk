rule bwa_se:
    input:
        fq = lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam = "{output_dir}/{sample}/{step}_bwa/bwa_{sample}.bam"
    params:
        index = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bwa/" +
            config[f"{wildcards.step}_ref_release"] + ".fa"
        ),
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bwa = config["bwa"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bwa/{sample}.log"
    threads: 10
    shell:
        """
        {params.bwa} mem -t {threads} {params.index} {input.fq} {params.other}  | \
        samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}

        total_reads=$(zcat {input.fq} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(100 * $aligned_reads / $total_reads, '.2f'))" )      
        echo "{wildcards.step} Alignment_Percentage $alignment_percentage%" >> {params.final_log}
        """

rule bwa_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bwa/bwa_{sample}.bam"
    params:
        index = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bwa/" +
            config[f"{wildcards.step}_ref_release"] + ".fa"
        ),
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bwa = config["bwa"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bwa/{sample}.log"
    threads: 10
    shell:
        """
        {params.bwa} mem -t {threads} {params.index} {input.fq1} {input.fq2}  {params.other}  | \
        samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}

        total_reads=$(zcat {input.fq1} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(50 * $aligned_reads / $total_reads, '.2f'))" )      
        echo "{wildcards.step} Alignment_Percentage $alignment_percentage%" >> {params.final_log}
        """

rule sp_bwa_se:
    input:
        fq = lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam = "{output_dir}/{sample}/{step}_bwa_spikein/bwa_spikein_{sample}.bam"
    params:
        index = config["ref_path"] + "/spikein/bwa/spikein.fa",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bwa = config["bwa"]
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bwa_spikein/{sample}.log"
    threads: 10
    shell:
        """
        {params.bwa} mem {params.other} -t {threads} {params.index} {input.fq} | \
        samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}
        """

rule sp_bwa_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bwa_spikein/bwa_spikein_{sample}.bam"
    params:
        index = config["ref_path"] + "/spikein/bwa/spikein.fa",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bwa = config["bwa"]
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bwa_spikein/{sample}.log"
    threads: 10
    shell:
        """
        {params.bwa} mem {params.other} -t {threads} {params.index} {input.fq1} {input.fq2} | \
        samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam}
        """


