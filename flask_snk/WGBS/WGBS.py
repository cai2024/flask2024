from os.path import join, exists

# 主要软件路径
fastp = "/data/biosoft/soft2024/fastp/fastp-0.23.4/fastp"
bowtie2="/data/biosoft/soft2024/bowtie2/bowtie2-2.5.2-linux-x86_64/bowtie2"
bowtie2_path="/data/biosoft/soft2024/bowtie2/bowtie2-2.5.2-linux-x86_64/"
bismark="/data/biosoft/soft2024/bismark/Bismark-0.24.2/bismark"
bismark_methylation_extractor="/data/biosoft/soft2024/bismark/Bismark-0.24.2/bismark_methylation_extractor"
picard="/data/biosoft/soft2024/picard/picard.jar"
qualimap="/data/biosoft/soft2024/qualimap/qualimap_v2.3/qualimap"
bedtools="/data/biosoft/soft2024/bedtools/bedtools2/bin/bedtools"
py_ref="/data/cailab/flask2024/flask_snk/ref"
java="/usr/bin/java"

spikein_ref="/data/reference2024/spikein/bismark"


# 加载配置文件
params_list = config['spikein_params'].split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)\t{param}_mc(CpG|CHG|CHH)' for param in params_list])


# 命令行加载
include: "trim.smk"
include: "mapping.smk"
include: "deduping.smk"
include: "call_modify.smk"
include: "cover.smk"
include: "logging.smk"



# 规则定义
rule all:
    input:
        expand("{output_dir}/bam/{align_soft}/{sample}_{mode}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"],align_soft=config['align_soft']),
        expand("{output_dir}/sp_bam/{align_soft}/sp_{sample}_{mode}.bam",
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"],align_soft=config['align_soft']),
        expand("{output_dir}/debam/{call_modify_soft}/{sample}_dedup_splitting_report.txt",
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"],call_modify_soft=config['call_modify_soft']),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        
              



