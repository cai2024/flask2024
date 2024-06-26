rule bowtie2_se: 
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft'])
    output:
        bam="{output_dir}/bam/bowtie2/{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{{sample}}_single.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} {params.other} ) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """

rule bowtie2_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/bam/bowtie2/{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{{sample}}_pair.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -1 {input.trim_fq1} -2 {input.trim_fq2} -S {params.sam}  -p {threads} {params.other} ) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """

rule sp_bowtie2_se:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft'])
    output:
        bam="{output_dir}/sp_bam/bowtie2/sp_{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bowtie2/spikein",
        sam="{output_dir}/sp_bam/sp_{{sample}}_single.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads}  {params.other}) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam} 
        samtools index {output.bam}
        rm {params.sam}
        """

rule sp_bowtie2_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/sp_bam/bowtie2/sp_{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bowtie2/spikein",
        sam="{output_dir}/sp_bam/sp_{{sample}}_pair.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -1 {input.trim_fq1} -2 {input.trim_fq2} -S {params.sam}  -p {threads} {params.other} ) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """





rule bwa_se:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft'])
    output:
        bam="{output_dir}/bam/bwa/{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/logs/bwa/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_genome"],
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.trim_fq1} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule bwa_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/bam/bwa/{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/logs/bwa/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_genome"],
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.trim_fq1} {input.trim_fq2}  2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule sp_bwa_se:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft'])
    output:
        bam="{output_dir}/sp_bam/bwa/sp_{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/logs/bwa/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bwa/spikein.fa",
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.trim_fq1} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule sp_bwa_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/sp_bam/bwa/sp_{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/logs/bwa/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bwa/spikein.fa",
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.trim_fq1} {input.trim_fq2} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

