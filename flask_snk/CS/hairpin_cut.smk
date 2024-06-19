rule hairpin_cut: 
    input:
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
    output:
        cut_fq="{output_dir}/trim/CGS_haipin/{{sample}}.fq".format(output_dir=config['output_dir']),
        cut_fq1="{output_dir}/trim/CGS_haipin/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
        cut_fq2="{output_dir}/trim/CGS_haipin/{{sample}}_cut_f2.fq".format(output_dir=config['output_dir']),                      
    log:
        "{output_dir}/logs/hairpin_cut/{{sample}}.log".format(output_dir=config['output_dir'])
    params:
        out_file="{output_dir}/trim/CGS_haipin/{{sample}}".format(output_dir=config['output_dir']),
        other=config['cut_params'],
    conda:
        "flask2024",
    threads: 20
    shell:
        """
        (python3 {py_ref}/new_hairpin_cut.py --fq1 {input.trim_fq1} --fq2 {input.trim_fq2}   \
        --outfile {params.out_file} {params.other} --parallel {threads}) > {log} 2>&1
        python3 {py_ref}/wild_stats.py --fq1 {output.cut_fq1} --fq2 {output.cut_fq2}
        """

