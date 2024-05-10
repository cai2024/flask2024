from os.path import join, exists

# 主要软件路径
fastp = "/data/biosoft/soft2024/fastp/fastp-0.23.4/fastp"
bowtie2="/data/biosoft/soft2024/bowtie2/bowtie2-2.5.2-linux-x86_64/bowtie2"
picard="/data/biosoft/soft2024/picard/picard.jar"
qualimap="/data/biosoft/soft2024/qualimap/qualimap_v2.3/qualimap"
py_ref="/data/cailab/flask_snk/ref"
# 加载配置文件
# 命令行加载

# 规则定义
rule all:
    input:
        expand("{output_dir}/bam/{sample}_{mode}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"]),

        expand("{output_dir}/debam/{sample}_dedup.bam",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/cover/{sample}/genome_results.txt",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),

rule fastp_se:
    input:
        fq1="{fq_in_path}/{{sample}}.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_single.fq.gz".format(output_dir=config['output_dir']) ,
        json="{output_dir}/trim/{{sample}}_single.json".format(output_dir=config['output_dir']),
        html="{output_dir}/trim/{{sample}}_single.html".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/fastp/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        config['fastp_params'],
    conda:
        "flask2024",
    threads: 2,
    shell:
        """
        {fastp} {params} -i {input.fq1} -o {output.trim_fq1} -j {output.json} -h {output.html} 2> {log}
        
        """

rule fastp_pe:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']) ,
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_{{sample}}_2.fq.gz".format(output_dir=config['output_dir']),
        json="{output_dir}/trim/{{sample}}_pair.json".format(output_dir=config['output_dir']),
        html="{output_dir}/trim/{{sample}}_pair.html".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/fastp/{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        config['fastp_params'],
    conda:
        "flask2024",  
    threads: 2,
    shell:
        """
        {fastp} {params} -i {input.fq1} -I {input.fq2} -o {output.trim_fq1} -O {output.trim_fq2} -j {output.json} -h {output.html} 2> {log}
        """
rule bowtie2_se: 
    input:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_single.fq.gz".format(output_dir=config['output_dir'])
    output:
        bam="{output_dir}/bam/{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{{sample}}_single.sam".format(output_dir=config['output_dir']),
        other=config['bowtie2_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """

rule bowtie2_pe:
    input:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_{{sample}}_2.fq.gz".format(output_dir=config['output_dir']),
    output:
        bam="{output_dir}/bam/{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{{sample}}_pair.sam".format(output_dir=config['output_dir']),
        other=config['bowtie2_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -1 {input.trim_fq1} -2 {input.trim_fq2} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """


rule mark_duplicate:
    input:
        bam="{output_dir}/bam/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"]),
    output:
        dedup_bam="{output_dir}/debam/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/duplication/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 4
    params:
        tmp_report="{output_dir}/debam/{{sample}}.txt".format(output_dir=config['output_dir']),
        other=config['picard_params']
    shell:
        """
        (/usr/bin/java -jar {picard} MarkDuplicatesWithMateCigar --REMOVE_DUPLICATES {params.other}-I {input.bam} -O {output.dedup_bam} -M {params.tmp_report}) 2> {log}
        samtools index {output.dedup_bam}        
        """

rule qualimap_bamqc:
    input: 
        dedup_bam="{output_dir}/debam/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    output: 
        "{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/{{sample}}".format(output_dir=config['output_dir']),
        java_mem_size="20G"
    threads: 10
    log:
        "{output_dir}/logs/qualimap/{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.dedup_bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size={params.java_mem_size}) > {log} 2>&1
        """

rule get_log:
    input:
        json="{output_dir}/trim/{{sample}}_{mode}.json".format(output_dir=config['output_dir'], mode=config["mode"]),
        cover="{output_dir}/cover/{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        bam="{output_dir}/bam/{{sample}}_{mode}.bam".format(output_dir=config['output_dir'],mode=config["mode"]),
        dedup_bam="{output_dir}/debam/{{sample}}_dedup.bam".format(output_dir=config['output_dir']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir'])
    conda:
        "flask2024",
    params:
        my_sample="{sample}"
    threads: 4

    shell:
        """
        name=$(echo {params.my_sample})
        raw_size=$(python3 {py_ref}/get_json.py {input.json} summary,before_filtering,total_bases)
        raw_read_num=$(python3 {py_ref}/get_json.py {input.json} summary,before_filtering,total_reads)
        dup_radio1=$(python3 {py_ref}/get_json.py {input.json} duplication,rate)

        clean_read_num=$(samtools view -c {input.bam}) 
        de_read_num=$(samtools view -c {input.dedup_bam})
        map_read_num=$(samtools view -c -F 4 {input.bam})
        
        cover=$(sed -n 's/.*There is a \([0-9.]*\)\\% of reference with a coverageData >= 1X.*/\\1/p' {input.cover})
        ref_size=$(sed -n 's/.*number of bases = \([0-9,]*\) bp.*/\\1/p' {input.cover} | tr -d ',')
        map_size=$(sed -n 's/.*number of mapped bases = \([0-9,]*\) bp.*/\\1/p' {input.cover} | tr -d ',')

        gc=$(python3 {py_ref}/get_json.py {input.json} summary,after_filtering,gc_content)
        q20=$(python3 {py_ref}/get_json.py {input.json} summary,after_filtering,q20_rate)
        q30=$(python3 {py_ref}/get_json.py {input.json} summary,after_filtering,q30_rate)
        chrm_num=$(samtools view -c {input.bam} chrM)



     
        [ "$(echo "$raw_read_num == 0" | bc)" -eq 1 ] && raw_read_num=100000000000
        [ "$(echo "$clean_read_num == 0" | bc)" -eq 1 ] && clean_read_num=100000000000
        [ "$(echo "$cover == 0" | bc)" -eq 1 ] && cover=100000000000
        [ "$(echo "$map_read_num == 0" | bc)" -eq 1 ] && map_read_num=100000000000
        ########################################
        sequence_size=$(echo "scale=4; $raw_size/1000000000" | bc)
        clean_read_ratio=$(echo "scale=4; $clean_read_num/$raw_read_num" | bc)
        duplication_rate=$(echo "scale=4; (($clean_read_num-$de_read_num)/$raw_read_num)+$dup_radio1" | bc)
        map_read_ratio=$(echo "scale=4; $map_read_num/$clean_read_num" | bc)
        average_depth=$(echo "scale=4; $map_size/($cover*$ref_size*0.01)" | bc)
        mit_ratio=$(echo "scale=4; $chrm_num/$map_read_num" | bc)

        gc_filter_ratio=$(echo "scale=4; $gc/1" | bc)
        q20_filter_ratio=$(echo "scale=4; $q20/1" | bc)
        q30_filter_ratio=$(echo "scale=4; $q30/1" | bc)


        echo -e "name\tsequence_size\tclean_read_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\tgc_filter_ratio\tmit_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "$name\t$sequence_size\t$clean_read_ratio\t$duplication_rate\t$map_read_ratio\t$cover\t$average_depth\t$gc_filter_ratio\t$mit_ratio\t$q20_filter_ratio\t$q30_filter_ratio\t$raw_read_num" >> {output.log}
        """






