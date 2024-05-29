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
spikein_ref = config['spikein_ref']



params_list = config['spikein_params'].split(',')
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
        expand("{output_dir}/bam/{call_modify_soft}/{sample}.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"],call_modify_soft=config["call_modify_soft"]),
        expand("{output_dir}/de_bam/{call_modify_soft}/{sample}_dedup.bam.all.bed.hmc",
               output_dir=config['output_dir'], sample=config["sampleList"],call_modify_soft=config["call_modify_soft"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
    














