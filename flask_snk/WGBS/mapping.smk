

rule bismark_se:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),       
    output:
        bam="{output_dir}/bam/bismark/{{sample}}_single.bam".format(output_dir=config['output_dir']),
        sort_bam="{output_dir}/bam/bismark/sort_{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bismark/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref_genome=config["ref_genome"],
        out_dir="{output_dir}/bam/bismark".format(output_dir=config['output_dir']),
        temp_bam="{output_dir}/bam/bismark/{{sample}}_trimmed_bismark_bt2.bam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        {bismark}  --path_to_bowtie2 {bowtie2_path} {params.ref_genome} {input.trim_fq1} {params.other} -o {params.out_dir} -p {threads}  > {log} 2>&1
        mv {params.temp_bam} {output.bam}
        samtools sort -@ {threads}  {output.bam} > {output.sort_bam}
        samtools index  {output.sort_bam}
        """

rule bismark_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/bam/bismark/{{sample}}_pair.bam".format(output_dir=config['output_dir']),
        sort_bam="{output_dir}/bam/bismark/sort_{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bismark/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref_genome=config["ref_genome"],
        out_dir="{output_dir}/bam/bismark".format(output_dir=config['output_dir']),
        temp_bam="{output_dir}/bam/bismark/{{sample}}_1_val_1_bismark_bt2_pe.bam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        {bismark} --path_to_bowtie2 {bowtie2_path} {params.ref_genome} -1 {input.trim_fq1} -2 {input.trim_fq2} {params.other} -o {params.out_dir} -p {threads}  > {log} 2>&1
        mv {params.temp_bam} {output.bam}
        samtools sort -@ {threads}  {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        """



rule sp_bismark_se:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'], trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/sp_bam/bismark/sp_{{sample}}_single.bam".format(output_dir=config['output_dir']),
        sort_bam="{output_dir}/sp_bam/bismark/sort_sp_{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/sp_bam/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bismark",
        out_dir="{output_dir}/sp_bam/bismark".format(output_dir=config['output_dir']),
        temp_bam="{output_dir}/sp_bam/bismark/{{sample}}_trimmed_bismark_bt2.bam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        {bismark} --path_to_bowtie2 {bowtie2_path} {params.ref} {input.trim_fq1} {params.other} -o {params.out_dir} -p {threads}  > {log} 2>&1
        mv {params.temp_bam} {output.bam}
        samtools sort -@ {threads} {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        """

rule sp_bismark_pe:
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/sp_bam/bismark/sp_{{sample}}_pair.bam".format(output_dir=config['output_dir']),
        sort_bam="{output_dir}/sp_bam/bismark/sort_sp_{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/sp_bam/sp_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=f"{ref_path}/spikein/bismark",
        out_dir="{output_dir}/sp_bam/bismark".format(output_dir=config['output_dir']),
        temp_bam="{output_dir}/sp_bam/bismark/{{sample}}_1_val_1_bismark_bt2_pe.bam".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        {bismark} --path_to_bowtie2 {bowtie2_path} {params.ref} -1 {input.trim_fq1} -2 {input.trim_fq2} {params.other} -o {params.out_dir} -p {threads}  > {log} 2>&1
        mv {params.temp_bam} {output.bam}
        samtools sort -@ {threads} {output.bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        """


