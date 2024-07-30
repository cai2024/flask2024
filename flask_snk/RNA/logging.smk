rule get_log:
    input: 
        bam="{output_dir}/{{sample}}/bam/{align_soft}/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"],align_soft=config['align_soft']),
        cover="{output_dir}/{{sample}}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        count_log="{output_dir}/{{sample}}/logs/count/{{sample}}.log".format(output_dir=config['output_dir']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir']),
    conda:
        "flask2024",
    params:
        raw_fq1="{fq_in_path}/{{sample}}".format(fq_in_path=config['fq_in_path']),
        trim_fq1="{output_dir}/{{sample}}/trim/{trim_soft}/{{sample}}".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        mode=config["mode"],
        my_sample="{sample}",
    threads: 4

    shell:
        """
        trim_fq1="{params.trim_fq1}$( [ {params.mode} = \'pair\' ] && echo \'_1_val_1\' || echo \'_trimmed\').fq.gz"
        size_clean=$(python3 {py_ref}/get_log.py size_clean --raw_fq1 {params.raw_fq1} --trim_fq1 {params.trim_fq1} --mode {params.mode})
        
        clean_num=$(seqtk seq $trim_fq1 | wc -l)
        map_num=$(samtools view -c -F 4 {input.bam})

        m1=$([ "{params.mode}" = "single" ] && echo "1" || echo "2" )
        map_ratio=$(python3 -c "print(format(4*int($map_num)/(int($clean_num)*int($m1)), '.4f'))")

        count=$(awk '/Successfully assigned alignments/ {{ print $6, $7 }}' {input.count_log})


        cover_depth=$(python3 {py_ref}/get_log.py cover_depth --input_file {input.cover} )
        
       

        mit_ratio=$(python3 {py_ref}/get_log.py mit_ratio --bam {input.bam} ) 
        gc_q20_q30=$(python3 {py_ref}/get_log.py gc_q20_q30 --fq_path {params.trim_fq1} --mode {params.mode} --bam {input.bam} --picard {picard} --java {java})
        raw_read_num=$(seqtk seq {params.raw_fq1}$( [ {params.mode} = \"pair\" ] && echo \"_1\" || echo \"\").fq.gz | wc -l | awk '{{print $1/4}}')
       




        echo -e "name\tsequence_size\tclean_read_ratio\tmap_read_ratio\tfeatureCounts_alignments\tcover_ratio\taverage_depth\tmit_ratio\tgc_filter_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "{params.my_sample}\t$size_clean\t$map_ratio\t$count\t$cover_depth\t$mit_ratio\t$gc_q20_q30\t$raw_read_num" >> {output.log}
        """

