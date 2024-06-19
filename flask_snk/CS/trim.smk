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
        other=config['trim_params'],
    conda:
        "flask2024",  
    threads: 3,
    shell:
        """
        {fastp} {params.other} -i {input.fq1} -I {input.fq2} -o {output.trim_fq1} -O {output.trim_fq2} -j {output.json} -h {output.html} -w {threads} > {log} 2>&1
        """


rule trim_galore_pe:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']),
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']),
    output:
        trim_fq1="{output_dir}/trim/trim_galore/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_galore/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/trim_galore/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        trim_dir="{output_dir}/trim/trim_galore/".format(output_dir=config['output_dir']),
        other=config['trim_params'],
    conda:
        "flask2024",
    threads: 3,
    shell:
        """
        trim_galore  {params.other} --paired --cores {threads} -o {params.trim_dir} {input.fq1} {input.fq2} > {log} 2>&1
        """


