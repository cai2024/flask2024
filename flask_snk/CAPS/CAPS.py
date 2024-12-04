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
bamUtil = config['bamUtil']



# 加载配置文件
params_list = config['spikein_params'].split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)' for param in params_list])

# 命令行加载

include: "../smk_modules/trim.smk"
include: "../smk_modules/mapping_bowtie2.smk"
include: "../smk_modules/mapping_bwa.smk"
include: "../smk_modules/deduping.smk"
include: "../smk_modules/de_cover.smk"

include: "logging.smk"
include: "../smk_modules/call_dup_astair.smk"
include: "../smk_modules/call_spikein_astair.smk"

# 规则定义
rule all:
    input:
        expand("{output_dir}/{sample}/bam/{align_soft}/{sample}_{mode}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"], align_soft=config['align_soft']),
        expand("{output_dir}/{sample}/de_bam/{dup_soft}/{sample}_dedup.bam",
               output_dir=config['output_dir'], sample=config["sampleList"], dup_soft=config['dup_soft']),


        expand("{output_dir}/{sample}/cover/{sample}/genome_results.txt",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/{sample}/sp_bam/{call_soft}/sp_{sample}_mCtoT_CpG.mods.gz",
               output_dir=config['output_dir'], sample=config["sampleList"],call_soft=config["call_soft"]),



        expand("{output_dir}/{sample}/sp_bam/{call_soft}/sp_{sample}_mCtoT_all.mods.gz",
               output_dir=config['output_dir'], sample=config["sampleList"],call_soft=config["call_soft"]),
        expand("{output_dir}/{sample}/de_bam/{call_soft}/de_{sample}_mCtoT_CpG.mods.gz",
               output_dir=config['output_dir'], sample=config["sampleList"],call_soft=config["call_soft"]),
        expand("{output_dir}/{sample}/de_bam/{call_soft}/de_{sample}_mCtoT_all.mods.gz",
               output_dir=config['output_dir'], sample=config["sampleList"],call_soft=config["call_soft"]),

        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
