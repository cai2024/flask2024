

rule bismark_call_methylation:
    input:
        bam="{output_dir}/debam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
    output:
        methylation_report="{output_dir}/debam/bismark/{{sample}}_dedup_splitting_report.txt".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/methylation/{{sample}}_methylation.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        temp_bam="{output_dir}/debam/{dup_soft}/temp_{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
        temp_report="{output_dir}/debam/bismark/temp_{{sample}}_dedup_splitting_report.txt".format(output_dir=config['output_dir']),
        out_dir="{output_dir}/debam/bismark".format(output_dir=config['output_dir']),
        call_modify_params=config['call_modify_params'],
    shell:
        """
        samtools sort -n {input.bam} > {params.temp_bam}
        {bismark_methylation_extractor} {params.temp_bam} {params.call_modify_params} --output {params.out_dir} > {log} 2>&1
        mv {params.temp_report} {output.methylation_report}
        """


rule bismark_call_methylation_sp:
    input:
        bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft'],mode=config['mode']),
    output:
        methylation_report="{output_dir}/sp_bam/bismark/sp_{{sample}}_{mode}.CX_report.txt".format(output_dir=config['output_dir'],mode=config['mode']),
    log:
        "{output_dir}/logs/sp_methylation/{{sample}}_methylation.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        out_dir="{output_dir}/sp_bam/bismark".format(output_dir=config['output_dir']),
        call_modify_params=config['call_modify_params'],
        ref=spikein_ref,
    shell:
        """
        {bismark_methylation_extractor} {input.bam} {params.call_modify_params} --output {params.out_dir}  --genome_folder {params.ref} --cytosine_report --merge_non_CpG --bedGraph --comprehensive --buffer_size 10G --CX --parallel {threads} > {log} 2>&1
        """

