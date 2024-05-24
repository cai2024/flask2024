rule get_log:
    input: 
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        bam="{output_dir}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        sort_bam="{output_dir}/bam/{align_soft}/sort_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        dedup_bam="{output_dir}/debam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
        sp_bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        methylation_report="{output_dir}/debam/{call_modify_soft}/{{sample}}_dedup_splitting_report.txt".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
        sp_mc="{output_dir}/sp_bam/{call_modify_soft}/sp_{{sample}}_{mode}.CX_report.txt".format(output_dir=config['output_dir'],mode=config['mode'],call_modify_soft=config['call_modify_soft']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir']),
    conda:
        "flask2024",
    params:
        fq_in_path=config['fq_in_path'],
        output_dir=config['output_dir'],
        trim_soft=config['trim_soft'],
        mode=config["mode"],
        my_sample="{sample}",
        spikein=config['spikein_params'],
    threads: 4

    shell:
        """
        size_clean=$(python3 {py_ref}/get_log.py size_clean --raw_fq1 {params.fq_in_path}/{params.my_sample} --trim_fq1 {params.output_dir}/trim/{params.trim_soft}/{params.my_sample} --mode {params.mode})
        dup_radio1=$(python3 {py_ref}/get_log.py get_json --json_file {params.output_dir}/trim/{params.trim_soft}/{params.my_sample}_{params.mode}.json --json_key duplication,rate)
        dup_map=$(python3 {py_ref}/get_log.py dup_map --bam {input.bam} --debam {input.dedup_bam} --dup_radio1 $dup_radio1)
        cover_depth=$(python3 {py_ref}/get_log.py cover_depth --input_file {input.cover} )
        all_spikein=$(python3 {py_ref}/get_log.py spikein_wgbs --bam {input.sp_bam} --spikein {params.spikein} --mc {input.sp_mc} --bedtools {bedtools})
        mit_ratio=$(python3 {py_ref}/get_log.py mit_ratio --bam {input.sort_bam} ) 
        gc_q20_q30=$(python3 {py_ref}/get_log.py gc_q20_q30 --fq_path {params.fq_in_path}/{params.my_sample} --mode {params.mode} --bam {input.bam} --picard {picard} --java {java})
        raw_read_num=$(seqtk seq {params.fq_in_path}/{params.my_sample}$( [ {params.mode} = \"pair\" ] && echo \"_1\" || echo \"\").fq.gz | wc -l | awk '{{print $1/4}}')
       

        mc_cpg=$(grep 'C methylated in CpG context' {input.methylation_report} | awk '{{print $NF}}')
        mc_chg=$(grep 'C methylated in CHG context' {input.methylation_report} | awk '{{print $NF}}')
        mc_chh=$(grep 'C methylated in CHH context' {input.methylation_report} | awk '{{print $NF}}')
        

        echo -e "name\tsequence_size\tclean_read_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\tmc_level(CpG|CHG|CHH)\t{spikein_list}\tmit_ratio\tgc_filter_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "{params.my_sample}\t$size_clean\t$dup_map\t$cover_depth\t$mc_cpg|$mc_chg|$mc_chh\t$all_spikein\t$mit_ratio\t$gc_q20_q30\t$raw_read_num" >> {output.log}
        """

