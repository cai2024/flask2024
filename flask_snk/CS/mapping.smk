rule bowtie2_se: 
    input:
        trim_fq1="{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir'])
    output:
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{align_soft}/{{sample}}.sam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """


rule bowtie2_se_sp:
    input:
        trim_fq1="{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir'])
    output:
        bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/bowtie2/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=spikein_ref,
        sam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.sam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        other=config['align_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """

