rule get_log:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']),
        trim_fq1="{output_dir}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        fq= "{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir']),
        cover="{output_dir}/cover/de_{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        dedup_bam="{output_dir}/debam/{align_soft}/de_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        hmc= "{output_dir}/debam/{align_soft}/de_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        mc="{output_dir}/debam/{align_soft}/de_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        sp_bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        sp_hmc= "{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        sp_mc="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
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
        name=$(echo {params.my_sample})
        size_clean=$(python3 {py_ref}/get_log.py size_clean --raw_fq1 {params.fq_in_path}/{params.my_sample} --trim_fq1 {params.output_dir}/trim/{params.trim_soft}/{params.my_sample} --mode {params.mode})
        match_model_ratio=$(echo "scale=4; $(cat {input.fq}  | wc -l)/$(zcat {input.trim_fq1} | wc -l)" | bc)
        dup_radio1=$(python3 {py_ref}/get_log.py get_json --json_file {params.output_dir}/trim/{params.trim_soft}/{params.my_sample}_{params.mode}.json --json_key duplication,rate)
        dup_map=$(python3 {py_ref}/get_log.py dup_map --bam {input.bam} --debam {input.dedup_bam} --dup_radio1 $dup_radio1)
        cover_depth=$(python3 {py_ref}/get_log.py cover_depth --input_file {input.cover} )
        all_spikein=$(python3 {py_ref}/get_log.py spikein --bam {input.sp_bam} --spikein {params.spikein} --hmc {input.sp_hmc} --mc {input.sp_mc} --bedtools {bedtools})
        mit_ratio=$(python3 {py_ref}/get_log.py mit_ratio --bam {input.bam} )
        gc_q20_q30=$(python3 {py_ref}/get_log.py gc_q20_q30 --fq_path {params.fq_in_path}/{params.my_sample} --mode {params.mode} --bam {input.bam} --picard {picard} --java {java})
        raw_read_num=$(seqtk seq {params.fq_in_path}/{params.my_sample}$( [ {params.mode} = \"pair\" ] && echo \"_1\" || echo \"\").fq.gz | wc -l | awk '{{print $1/4}}')
        
        
        hmc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc})
        mc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})
        mc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc}) 
        mc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})

        hmc="0"$(echo "scale=4; $hmc_cpg/1" | bc)"|0"$(echo "scale=4; $hmc_chg/1" | bc)"|0"$(echo "scale=4; $hmc_chh/1" | bc)
        mc="0"$(echo "scale=4; $mc_cpg/1" | bc)"|0"$(echo "scale=4; $mc_chg/1" | bc)"|0"$(echo "scale=4; $mc_chh/1" | bc)
        




     


        echo -e "name\tsequence_size\tclean_read_ratio\tmatch_model_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\thmc_level(CpG|CHG|CHH)\tmc_level(CpG|CHG|CHH)\t{spikein_list}\tmit_ratio\tgc_filter_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "$name\t$size_clean\t$match_model_ratio\t$dup_map\t$cover_depth\t$hmc\t$mc\t$all_spikein\t$mit_ratio\t$gc_q20_q30\t$raw_read_num" >> {output.log}
        """




