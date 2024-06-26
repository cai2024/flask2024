

rule bismark_call_methylation:
    input:
        bam="{output_dir}/de_bam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
    output:
        methylation_report="{output_dir}/de_bam/bismark/{{sample}}_dedup_splitting_report.txt".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/methylation/{{sample}}_methylation.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        out_dir="{output_dir}/de_bam/bismark".format(output_dir=config['output_dir']),
        other=config['call_modify_params'],
    shell:
        """
        {bismark_methylation_extractor} {input.bam} {params.other} --output {params.out_dir} > {log} 2>&1
        """


rule sp_bismark_call_methylation:
    input:
        bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft'],mode=config['mode']),
    output:
        methylation_report="{output_dir}/sp_bam/bismark/sp_{{sample}}_{mode}.CX_report.txt".format(output_dir=config['output_dir'],mode=config['mode']),
    log:
        "{output_dir}/logs/sp_methylation/sp_{{sample}}_methylation.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        out_dir="{output_dir}/sp_bam/bismark".format(output_dir=config['output_dir']),
        ref=f"{ref_path}/spikein",
        other=config['call_modify_params'],
    shell:
        """
        {bismark_methylation_extractor} {input.bam} {params.other} --output {params.out_dir}  --genome_folder {params.ref} --cytosine_report --merge_non_CpG --bedGraph --comprehensive --buffer_size 10G --CX --parallel {threads} > {log} 2>&1
        """

