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




