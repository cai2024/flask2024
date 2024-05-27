rule get_log:
    input:
        cut_fq="{output_dir}/trim/{cut_soft}/{{sample}}.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        de_bam="{output_dir}/de_bam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
        sp_bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        sp_hmc= "{output_dir}/sp_bam/{call_modify_soft}/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
        sp_mc="{output_dir}/sp_bam/{call_modify_soft}/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        hmc= "{output_dir}/de_bam/{call_modify_soft}/{{sample}}_dedup.bam.all.bed.hmc".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
        mc="{output_dir}/de_bam/{call_modify_soft}/{{sample}}_dedup.bam.all.bed.mc".format(output_dir=config['output_dir'],call_modify_soft=config['call_modify_soft']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir']),
    conda:
        "flask2024",
    params:
        raw_fq1="{fq_in_path}/{{sample}}".format(fq_in_path=config['fq_in_path']),
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        json="{output_dir}/trim/{trim_soft}/{{sample}}_{mode}.json".format(output_dir=config['output_dir'],trim_soft=config['trim_soft'],mode=config["mode"]),
        spikein=config['spikein_params'],
        my_sample="{sample}",
    threads: 4
    shell:
        """
        size_clean=$(python3 {py_ref}/get_log.py size_clean --raw_fq1 {params.raw_fq1} --trim_fq1 {params.trim_fq1} --mode pair)
        match_model_ratio=$(echo "scale=4; $(cat {input.cut_fq}  | wc -l)/$(zcat {params.trim_fq1}_1_val_1.fq.gz  | wc -l)" | bc)
        dup_radio1=$(python3 {py_ref}/get_log.py get_json --json_file {params.json} --json_key duplication,rate)
        dup_map=$(python3 {py_ref}/get_log.py dup_map --bam {input.bam} --debam {input.de_bam} --dup_radio1 $dup_radio1)
        cover_depth=$(python3 {py_ref}/get_log.py cover_depth --input_file {input.cover} )
        all_spikein=$(python3 {py_ref}/get_log.py spikein --bam {input.sp_bam} --spikein {params.spikein} --hmc {input.sp_hmc} --mc {input.sp_mc} --bedtools {bedtools})
        mit_ratio=$(python3 {py_ref}/get_log.py mit_ratio --bam {input.bam} ) 
        gc_q20_q30=$(python3 {py_ref}/get_log.py gc_q20_q30 --fq_path {params.trim_fq1} --mode pair --bam {input.bam} --picard {picard} --java {java})
        raw_read_num=$(seqtk seq {params.raw_fq1}_1.fq.gz | wc -l | awk '{{print $1/4}}')
  
        hmc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc})
        mc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})
        mc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc}) 
        mc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})

        hmc="0"$(echo "scale=4; $hmc_cpg/1" | bc)"|0"$(echo "scale=4; $hmc_chg/1" | bc)"|0"$(echo "scale=4; $hmc_chh/1" | bc)
        mc="0"$(echo "scale=4; $mc_cpg/1" | bc)"|0"$(echo "scale=4; $mc_chg/1" | bc)"|0"$(echo "scale=4; $mc_chh/1" | bc)
        

        echo -e "name\tsequence_size\tclean_read_ratio\tmatch_model_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\thmc_level(CpG|CHG|CHH)\tmc_level(CpG|CHG|CHH)\t{spikein_list}\tmit_ratio\tgc_filter_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "{params.my_sample}\t$size_clean\t$match_model_ratio\t$dup_map\t$cover_depth\t$hmc\t$mc\t$all_spikein\t$mit_ratio\t$gc_q20_q30\t$raw_read_num" >> {output.log}
        """




