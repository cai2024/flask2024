from os.path import join, exists

# 主要软件路径
fastp = "/data/biosoft/soft2024/fastp/fastp-0.23.4/fastp"
bowtie2="/data/biosoft/soft2024/bowtie2/bowtie2-2.5.2-linux-x86_64/bowtie2"
picard="/data/biosoft/soft2024/picard/picard.jar"
qualimap="/data/biosoft/soft2024/qualimap/qualimap_v2.3/qualimap"
bedtools="/data/biosoft/soft2024/bedtools/bedtools2/bin/bedtools"
samtools="/data/biosoft/software/samtools-1.9/samtools"
# 加载配置文件
py_ref="/data/cailab/flask_snk/ref"
spikein_ref="/data/reference2024/spikein/bowtie2/spikein"
log_params = config['spikein_params']
params_list = log_params.split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)\t{param}_hmc(CpG|CHG|CHH)\t{param}_mc(CpG|CHG|CHH)' for param in params_list])
# 命令行加载

# 规则定义
rule all:
    input:
        expand("{output_dir}/bam/{sample}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/bam/{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/debam/de_{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/cover/de_{sample}/genome_results.txt",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/sp_bam/sp_{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"]),

rule fastp_pe:
    input:
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']) ,
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_{{sample}}_2.fq.gz".format(output_dir=config['output_dir']),
        json="{output_dir}/trim/{{sample}}.json".format(output_dir=config['output_dir']),
        html="{output_dir}/trim/{{sample}}.html".format(output_dir=config['output_dir']),
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

rule hairpin_cut: 
    input:
        trim_fq1="{output_dir}/trim/trim_{{sample}}_1.fq.gz".format(output_dir=config['output_dir']),
        trim_fq2="{output_dir}/trim/trim_{{sample}}_2.fq.gz".format(output_dir=config['output_dir']),

    output:
        fq="{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
        fq2="{output_dir}/trim/{{sample}}_cut_f2.fq".format(output_dir=config['output_dir']),
                      
    log:
        "{output_dir}/logs/hairpin_cut/{{sample}}.log".format(output_dir=config['output_dir'])
    params:
        out_file="{output_dir}/trim/{{sample}}".format(output_dir=config['output_dir']),
        other=config['hairpin_cut_params'],
    conda:
        "flask2024",
    threads: 20
    shell:
        """
        (python3 {py_ref}/new_hairpin_cut.py --fq1 {input.trim_fq1} --fq2 {input.trim_fq2}   \
        --outfile {params.out_file} {params.other} --parallel {threads}) > {log} 2>&1
        python3 {py_ref}/wild_stats.py --fq1 {output.fq1} --fq2 {output.fq2}
        """




rule bowtie2_se: 
    input:
        trim_fq1="{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir'])
    output:
        bam="{output_dir}/bam/{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=config["ref_genome"],
        sam="{output_dir}/bam/{{sample}}.sam".format(output_dir=config['output_dir']),
        other=config['bowtie2_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """


rule extract_mc_hmc:
    input:
        bam="{output_dir}/bam/{{sample}}.bam".format(output_dir=config['output_dir']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/bam/{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/bam/{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=config["ref_genome"],
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """


rule mark_duplicate: 
    input:
        bam="{output_dir}/bam/{{sample}}.bam".format(output_dir=config['output_dir']),
        fq1="{fq_in_path}/{{sample}}_1.fq.gz".format(fq_in_path=config['fq_in_path']) ,
        fq2="{fq_in_path}/{{sample}}_2.fq.gz".format(fq_in_path=config['fq_in_path']) ,
    output:
        temp_bam=temp("{output_dir}/debam/de_temp_{{sample}}.bam".format(output_dir=config['output_dir'])),
        bam="{output_dir}/debam/de_{{sample}}.bam".format(output_dir=config['output_dir'])
    log:
        "{output_dir}/logs/dedup/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 10
    shell:
        """
        (python3 {py_ref}/remove_dup.py --input {input.bam} --output {output.temp_bam} --fq1 {input.fq1} --fq2 {input.fq2} ) 2> {log}
        samtools sort {output.temp_bam} > {output.bam}
        samtools index {output.bam}
        """

rule de_extract_mc_hmc:
    input:
        bam="{output_dir}/debam/de_{{sample}}.bam".format(output_dir=config['output_dir']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/debam/de_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/debam/de_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/de_extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=config["ref_genome"]
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """




rule bowtie2_se_sp:
    input:
        trim_fq1="{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir'])
    output:
        bam="{output_dir}/sp_bam/sp_{{sample}}.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/bowtie2/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    params:
        ref=spikein_ref,
        sam="{output_dir}/sp_bam/sp_{{sample}}.sam".format(output_dir=config['output_dir']),
        other=config['bowtie2_params'],
    shell:
        """
        ({bowtie2} -x {params.ref} -U {input.trim_fq1} -S {params.sam}  -p {threads} ) 2> {log} 
        samtools view -b {params.sam} | samtools sort > {output.bam}
        samtools index {output.bam}
        rm {params.sam}
        """


rule extract_mc_hmc_sp:
    input:
        bam="{output_dir}/sp_bam/sp_{{sample}}.bam".format(output_dir=config['output_dir']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/sp_bam/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/sp_bam/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=spikein_ref
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """


















rule qualimap_bamqc:
    input: 
        dedup_bam="{output_dir}/debam/de_{{sample}}.bam".format(output_dir=config['output_dir']),
    output: 
        "{output_dir}/cover/de_{{sample}}/genome_results.txt".format(output_dir=config['output_dir'])
    params:
        out_dir="{output_dir}/cover/de_{{sample}}".format(output_dir=config['output_dir']),
        java_mem_size="20G"
    threads: 10
    log:
        "{output_dir}/logs/qualimap/de_{{sample}}.log".format(output_dir=config['output_dir']),
    shell:
        """
        ({qualimap} bamqc -bam {input.dedup_bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size={params.java_mem_size}) > {log} 2>&1
        """




rule get_log:
    input:
        fq= "{output_dir}/trim/{{sample}}.fq".format(output_dir=config['output_dir']),
        json="{output_dir}/trim/{{sample}}.json".format(output_dir=config['output_dir'], mode=config["mode"]),
        cover="{output_dir}/cover/de_{{sample}}/genome_results.txt".format(output_dir=config['output_dir']),
        bam="{output_dir}/bam/{{sample}}.bam".format(output_dir=config['output_dir'],mode=config["mode"]),
        dedup_bam="{output_dir}/debam/de_{{sample}}.bam".format(output_dir=config['output_dir']),
        hmc= "{output_dir}/debam/de_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/debam/de_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
        sp_bam="{output_dir}/sp_bam/sp_{{sample}}.bam".format(output_dir=config['output_dir']),
        sp_hmc= "{output_dir}/sp_bam/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        sp_mc="{output_dir}/sp_bam/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    output:
        log="{output_dir}/final_log/{{sample}}.log".format(output_dir=config['output_dir']),
        bed="{output_dir}/sp_bam/sp_{{sample}}.bed".format(output_dir=config['output_dir']),
    conda:
        "flask2024",
    params:
        my_sample="{sample}",
        spikein=config['spikein_params'],
    threads: 4

    shell:
        """
        name=$(echo {params.my_sample})
        raw_size=$(python3 {py_ref}/get_json.py {input.json} summary,before_filtering,total_bases)
        raw_read_num=$(python3 {py_ref}/get_json.py {input.json} summary,before_filtering,total_reads)
        trimed_read_num=$(python3 {py_ref}/get_json.py {input.json} summary,after_filtering,total_reads)
        model_read_num=$(wc -l < {input.fq})
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
        
        hmc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc}) 
        hmc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.hmc})
        mc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})
        mc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc}) 
        mc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {input.mc})

        hmc="0"$(echo "scale=4; $hmc_cpg/1" | bc)"|0"$(echo "scale=4; $hmc_chg/1" | bc)"|0"$(echo "scale=4; $hmc_chh/1" | bc)
        mc="0"$(echo "scale=4; $mc_cpg/1" | bc)"|0"$(echo "scale=4; $mc_chg/1" | bc)"|0"$(echo "scale=4; $mc_chh/1" | bc)
        {bedtools} genomecov -ibam {input.sp_bam}  -bga > {output.bed}
        all_spikein=$(python3 {py_ref}/spikein.py --bed {output.bed} --spikein {params.spikein}  --hmc {input.sp_hmc} --mc {input.sp_mc})
        


     
        [ "$(echo "$raw_read_num == 0" | bc)" -eq 1 ] && raw_read_num=100000000000
        [ "$(echo "$clean_read_num == 0" | bc)" -eq 1 ] && clean_read_num=100000000000
        [ "$(echo "$cover == 0" | bc)" -eq 1 ] && cover=100000000000
        [ "$(echo "$map_read_num == 0" | bc)" -eq 1 ] && map_read_num=100000000000
        ########################################
        sequence_size=$(echo "scale=4; $raw_size/1000000000" | bc)
        clean_read_ratio=$(echo "scale=4; $trimed_read_num/$raw_read_num" | bc)
        match_model_ratio=$(echo "scale=4; $model_read_num/($trimed_read_num*2)" | bc)
        duplication_rate=$(echo "scale=4; (($clean_read_num-$de_read_num)/$clean_read_num)" | bc)
        map_read_ratio=$(echo "scale=4; $map_read_num/$clean_read_num" | bc)
        average_depth=$(echo "scale=4; $map_size/($cover*$ref_size*0.01)" | bc)


        mit_ratio=$(echo "scale=4; $chrm_num/$map_read_num" | bc)
        gc_filter_ratio=$(echo "scale=4; $gc/1" | bc)
        q20_filter_ratio=$(echo "scale=4; $q20/1" | bc)
        q30_filter_ratio=$(echo "scale=4; $q30/1" | bc)


        echo -e "name\tsequence_size\tclean_read_ratio\tmatch_model_ratio\tduplication_rate\tmap_read_ratio\tcover_ratio\taverage_depth\thmc_level(CpG|CHG|CHH)\tmc_level(CpG|CHG|CHH)\t{spikein_list}\tgc_filter_ratio\tmit_ratio\tq20_filter_ratio\tq30_filter_ratio\ttotal_reads_raw" > {output.log}

        echo -e "$name\t$sequence_size(G)\t$clean_read_ratio\t$match_model_ratio\t$duplication_rate\t$map_read_ratio\t$cover\t$average_depth\t$hmc\t$mc\t$all_spikein\t$gc_filter_ratio\t$mit_ratio\t$q20_filter_ratio\t$q30_filter_ratio\t$raw_read_num" >> {output.log}
        """








