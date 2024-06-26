from os.path import join, exists



fastp = config['fastp']
bowtie2 = config['bowtie2']
bowtie2_path = config['bowtie2_path']
bwa = config['bwa']
bismark = config['bismark']
bismark_methylation_extractor = config['bismark_methylation_extractor']
picard = config['picard']
qualimap = config['qualimap']
bedtools = config['bedtools']
py_ref = config['py_ref']
java = config['java']
ref_path = config['ref_path']




# 加载配置文件
params_list = config['spikein_params'].split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)' for param in params_list])

# 命令行加载

include: "trim.smk"
include: "mapping.smk"
include: "deduping.smk"
include: "cover.smk"
include: "logging.smk"



# 规则定义
rule all:
    input:
        expand("{output_dir}/bam/{align_soft}/{sample}_{mode}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"], align_soft=config['align_soft']),
        expand("{output_dir}/de_bam/{dup_soft}/{sample}_dedup.bam",
               output_dir=config['output_dir'], sample=config["sampleList"], dup_soft=config['dup_soft']),


        expand("{output_dir}/cover/{sample}/genome_results.txt",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),

