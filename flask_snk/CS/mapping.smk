rule bowtie2_se: 
    input:
        cut_fq="{output_dir}/{{sample}}/trim/{cut_soft}/{{sample}}.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        bam="{output_dir}/{{sample}}/bam/bowtie2/{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/{{sample}}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] +"/bowtie2/" + config["ref_release"],
        sam="{output_dir}/{{sample}}/bam/{{sample}}.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.cut_fq} -S {params.sam}  -p {threads} {params.other} ) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """


rule sp_bowtie2_se:
    input:
        cut_fq="{output_dir}/{{sample}}/trim/{cut_soft}/{{sample}}.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        bam="{output_dir}/{{sample}}/sp_bam/bowtie2/sp_{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/{{sample}}/logs/bowtie2/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bowtie2/spikein",
        sam="{output_dir}/{{sample}}/sp_bam/sp_{{sample}}.sam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.cut_fq} -S {params.sam}  -p {threads}  {params.other}) > {log} 2>&1
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam} 
        samtools index {output.bam}
        rm {params.sam}
        """




rule bwa_se:
    input:
        cut_fq="{output_dir}/{{sample}}/trim/{cut_soft}/{{sample}}.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        bam="{output_dir}/{{sample}}/bam/bwa/{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/{{sample}}/logs/bwa/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] +"/bowtie2/" + config["ref_release"] +".fa",
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.cut_fq} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """



rule sp_bwa_se:
    input:
        cut_fq="{output_dir}/{{sample}}/trim/{cut_soft}/{{sample}}.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        bam="{output_dir}/{{sample}}/sp_bam/bwa/sp_{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        log="{output_dir}/{{sample}}/logs/bwa/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bwa/spikein.fa",
        other=config['align_params'],
    shell:
        """
        {bwa} mem {params.other} -t {threads} {params.ref} {input.cut_fq} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


