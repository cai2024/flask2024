rule bismark_se:
    input:
        fq = lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam = "{output_dir}/{sample}/{step}_bismark/bismark_{sample}.bam",
    params:
        index = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bismark/"
        ),
        out_dir = "{output_dir}/{sample}/{step}_bismark/",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bismark = config["bismark"],
        bowtie2_path = config["bowtie2_path"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bismark/{sample}.log"
    threads: 10
    shell:
        """
        {params.bismark} --path_to_bowtie2 {params.bowtie2_path} {params.index} {input.fq} {params.other} -o {params.out_dir} -p {threads} > {log} 2>&1 
        if grep -q "ERR" {log}; then echo "Bismark encountered an error，high likelihood of memory error!" >&2; exit 1; fi 
        temp_bam=$(basename {input.fq} .fq.gz)_bismark_bt2.bam
        samtools view -bh -F 256 -F 2048 -@ {threads} {params.out_dir}$temp_bam |  samtools sort -@ {threads}  > {output.bam}
        samtools index {output.bam}



        total_reads=$(zcat {input.fq} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(100 * $aligned_reads / $total_reads, '.2f'))" )
        echo "{wildcards.step} Alignment_Percentage $alignment_percentage%" >> {params.final_log}

        """

rule bismark_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bismark/bismark_{sample}.bam",
    params:
        index = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/bismark/"
        ),
        out_dir = "{output_dir}/{sample}/{step}_bismark/",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bismark = config["bismark"],
        bowtie2_path = config["bowtie2_path"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bismark/{sample}.log"
    threads: 10
    shell:
        """
        {params.bismark} --path_to_bowtie2 {params.bowtie2_path} {params.index} -1 {input.fq1} -2 {input.fq2} {params.other} -o {params.out_dir} -p {threads} > {log} 2>&1 \
        && if grep -q "ERR" {log}; then echo "Bismark encountered an error, high likelihood of memory error!" >&2; exit 1; fi \
        && temp_bam=$(basename {input.fq1} .fq.gz)_bismark_bt2_pe.bam \
        && samtools view -bh -F 256 -F 2048 -@ {threads} {params.out_dir}$temp_bam | samtools sort -@ {threads}  > {output.bam} \
        && samtools index {output.bam}

        total_reads=$(zcat {input.fq1} | echo $((`wc -l` / 4)))
        aligned_reads=$(samtools view -c -F 4 {output.bam})
        alignment_percentage=$(python3 -c "print(format(50 * $aligned_reads / $total_reads, '.2f'))" )      
        echo "{wildcards.step} Alignment_Percentage $alignment_percentage%" >> {params.final_log}
        """

rule sp_bismark_se:
    input:
        fq = lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        bam = "{output_dir}/{sample}/{step}_bismark_spikein/bismark_spikein_{sample}.bam",
    params:
        index = config["ref_path"] + "/spikein/bismark/",
        out_dir = "{output_dir}/{sample}/{step}_bismark_spikein/",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bismark = config["bismark"],
        bowtie2_path = config["bowtie2_path"],
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bismark_spikein/{sample}.log"
    threads: 10
    shell:
        """
        {params.bismark} --path_to_bowtie2 {params.bowtie2_path} {params.index} {input.fq} {params.other} -o {params.out_dir} -p {threads} > {log} 2>&1 \
        && if grep -q "ERR" {log}; then echo "Bismark encountered an error，high likelihood of memory error!" >&2; exit 1; fi \
        && temp_bam=$(basename {input.fq} .fq.gz)_bismark_bt2.bam \
        && samtools sort -@ {threads} {params.out_dir}$temp_bam > {output.bam} \
        && samtools index {output.bam}
        """

rule sp_bismark_pe:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        bam = "{output_dir}/{sample}/{step}_bismark_spikein/bismark_spikein_{sample}.bam",
    params:
        index = config["ref_path"] + "/spikein/bismark/",
        out_dir = "{output_dir}/{sample}/{step}_bismark_spikein/",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        bismark = config["bismark"],
        bowtie2_path = config["bowtie2_path"],
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_bismark_spikein/{sample}.log"
    threads: 10
    shell:
        """
        {params.bismark} --path_to_bowtie2 {params.bowtie2_path} {params.index} -1 {input.fq1} -2 {input.fq2} {params.other} -o {params.out_dir} -p {threads} > {log} 2>&1 \
        && if grep -q "ERR" {log}; then echo "Bismark encountered an error, high likelihood of memory error!" >&2; exit 1; fi \
        && temp_bam=$(basename {input.fq1} .fq.gz)_bismark_bt2_pe.bam \
        && samtools sort -@ {threads} {params.out_dir}$temp_bam > {output.bam} \
        && samtools index {output.bam}
        """
