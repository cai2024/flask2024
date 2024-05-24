rule fastp_se:
    input:
        fq1="{fq_in_path}/{{sample}}.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        trim_fq1="{output_dir}/trim/fastp/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir']) ,
        json="{output_dir}/trim/fastp/{{sample}}_single.json".format(output_dir=config['output_dir']),
        html="{output_dir}/trim/fastp/{{sample}}_single.html".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/fastp/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        config['trim_params'],
    conda:
        "flask2024",
    threads: 3,
    shell:
        """
        {fastp} {params} -i {input.fq1} -o {output.trim_fq1} -j {output.json} -h {output.html} -w {threads}  2> {log}
        """

rule fastp_pe:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']) ,
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        trim_fq1="{output_dir}/trim/fastp/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/fastp/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir']),
        json="{output_dir}/trim/fastp/{{sample}}_pair.json".format(output_dir=config['output_dir']),
        html="{output_dir}/trim/fastp/{{sample}}_pair.html".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/fastp/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        config['trim_params'],
    conda:
        "flask2024",  
    threads: 3,
    shell:
        """
        {fastp} {params} -i {input.fq1} -I {input.fq2} -o {output.trim_fq1} -O {output.trim_fq2} -j {output.json} -h {output.html} -w {threads} 2> {log}
        """


rule trim_galore_se:
    input:
        fq1="{fq_in_path}/{{sample}}.fq.gz".format(fq_in_path=config['fq_in_path']),
    output:
        trim_fq="{output_dir}/trim/trim_galore/trim_{{sample}}.fq.gz".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/trim_galore/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        extra=config['trim_params'],
    conda:
        "flask2024",  
    threads: 3,
    shell:
        """
        trim_galore {params.extra} --cores {threads} -o {output.trim_fq|dirname} {input.fq1} > {output.report} 2> {log}
        """
rule trim_galore_pe:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']),
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']),
    output:
        trim_fq1="{output_dir}/trim/trim_galore/trim_{{sample}}_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_galore/trim_{{sample}}_2.fq.gz".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/trim_galore/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        extra=config['trim_params'],
    conda:
        "flask2024",
    threads: 3,
    shell:
        """
        trim_galore {params.extra} --paired --cores {threads} -o {output.trim_fq1|dirname} {input.fq1} {input.fq2} > {output.report} 2> {log}
        """


