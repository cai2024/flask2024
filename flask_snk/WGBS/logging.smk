rule get_log:
    input: 
        sort_bam="{output_dir}/bam/{align_soft}/sort_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        de_sort_bam="{output_dir}/de_bam/{dup_soft}/sort_{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
        sp_sort_bam="{output_dir}/sp_bam/{align_soft}/sort_sp_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        sp_bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        mc_report="{output_dir}/de_bam/{call_modify_soft}/{{sample}}_dedup_splitting_report.txt".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
        mc_sp="{output_dir}/sp_bam/{call_modify_soft}/sp_{{sample}}_{mode}.CX_report.txt".format(output_dir=config['output_dir'],mode=config['mode'],call_modify_soft=config['call_modify_soft']),
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir']),
    conda:
        "flask2024",
    params: 
        raw_fq1="{fq_in_path}/{{sample}}".format(fq_in_path=config['fq_in_path']),
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        json="{output_dir}/trim/{trim_soft}/{{sample}}_{mode}.json".format(output_dir=config['output_dir'],trim_soft=config['trim_soft'],mode=config["mode"]),

        spikein=config['spikein_params'],
        mode=config["mode"],
        my_sample="{sample}",

    threads: 4
    shell:
        """
        size_clean=$(python3 {py_ref}/get_log.py size_clean --raw_fq1 {params.raw_fq1} --trim_fq1 {params.trim_fq1} --mode {params.mode})
        dup_radio1=$(python3 {py_ref}/get_log.py get_json --json_file {params.json} --json_key duplication,rate)
        dup_map=$(python3 {py_ref}/get_log.py dup_map --bam {input.sort_bam} --debam {input.de_sort_bam} --dup_radio1 $dup_radio1)
        cover_depth=$(python3 {py_ref}/get_log.py cover_depth --input_file {input.cover} )
        all_spikein=$(python3 {py_ref}/get_log.py spikein_wgbs --bam {input.sp_bam} --sort_bam {input.sp_sort_bam} --mc {input.mc_sp} --spikein {params.spikein} --bedtools {bedtools})
        mit_ratio=$(python3 {py_ref}/get_log.py mit_ratio --bam {input.sort_bam} ) 
        gc_q20_q30=$(python3 {py_ref}/get_log.py gc_q20_q30 --fq_path {params.trim_fq1} --mode {params.mode} --bam {input.sort_bam} --picard {picard} --java {java})
        raw_read_num=$(seqtk seq {params.raw_fq1}$( [ {params.mode} = \"pair\" ] && echo \"_1\" || echo \"\").fq.gz | wc -l | awk '{{print $1/4}}')


       

        mc_cpg=$(grep 'C methylated in CpG context' {input.mc_report} | awk '{{print $NF}}')
        mc_chg=$(grep 'C methylated in CHG context' {input.mc_report} | awk '{{print $NF}}')
        mc_chh=$(grep 'C methylated in CHH context' {input.mc_report} | awk '{{print $NF}}')
        

        echo -e "name\tsequence_size\tclean_read_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\tmc_level(CpG|CHG|CHH)\t{spikein_list}\tmit_ratio\tgc_filter_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "{params.my_sample}\t$size_clean\t$dup_map\t$cover_depth\t$mc_cpg|$mc_chg|$mc_chh\t$all_spikein\t$mit_ratio\t$gc_q20_q30\t$raw_read_num" >> {output.log}
        """

