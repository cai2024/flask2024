from os.path import join, exists

# 主要软件路径
fastp = "/data/biosoft/soft2024/fastp/fastp-0.23.4/fastp"
bowtie2="/data/biosoft/soft2024/bowtie2/bowtie2-2.5.2-linux-x86_64/bowtie2"
picard="/data/biosoft/soft2024/picard/picard.jar"
qualimap="/data/biosoft/soft2024/qualimap/qualimap_v2.3/qualimap"
bedtools="/data/biosoft/soft2024/bedtools/bedtools2/bin/bedtools"
samtools="/data/biosoft/software/samtools-1.9/samtools"
java="/usr/bin/java"
# 加载配置文件
py_ref="/data/cailab/flask2024/flask_snk/ref"
spikein_ref="/data/reference2024/spikein/bowtie2/spikein"
log_params = config['spikein_params']
params_list = log_params.split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)\t{param}_hmc(CpG|CHG|CHH)\t{param}_mc(CpG|CHG|CHH)' for param in params_list])
# 命令行加载
include: "trim.smk"
include: "mapping.smk"
include: "deduping.smk"
include: "cover.smk"
include: "logging.smk"
include: "hairpin_cut.smk"
include: "extract_mc_hmc.smk"


# 规则定义
rule all:
    input:
        expand("{output_dir}/bam/{align_soft}/{sample}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"],align_soft=config['align_soft']),
        expand("{output_dir}/bam/{align_soft}/{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"],align_soft=config["align_soft"]),
        expand("{output_dir}/debam/{align_soft}/de_{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"],align_soft=config["align_soft"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
    














