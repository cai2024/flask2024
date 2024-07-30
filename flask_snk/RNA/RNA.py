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
STAR = config['STAR']
featureCounts  = config['featureCounts']
ref_path = config['ref_path']
ref_species = config["ref_species"]
ref_release = config["ref_release"]
gtf_path = ref_path + "/" + ref_species + "/" + ref_release + "/" + ref_release + ".gtf" 



# 加载配置文件
params_list = config['spikein_params'].split(',')
spikein_list = '\t'.join([f'{param}_(cover|dep)' for param in params_list])

# 命令行加载

include: "../smk_modules/trim.smk"
include: "../smk_modules/mapping_star.smk"
include: "../smk_modules/mapping_bowtie2.smk"
include: "../smk_modules/mapping_bwa.smk"
include: "../smk_modules/cover.smk"
include: "../smk_modules/count_featureCounts.smk"
include: "logging.smk"



# 规则定义
rule all:
    input:
        expand("{output_dir}/{sample}/bam/{align_soft}/{sample}_{mode}.bam", 
               output_dir=config['output_dir'], sample=config["sampleList"], mode=config["mode"], align_soft=config['align_soft']),
 
        expand("{output_dir}/{sample}/cover/{sample}/genome_results.txt",
               output_dir=config['output_dir'], sample=config["sampleList"]),
        expand("{output_dir}/{sample}/count/{count_soft}/{sample}_expression_matrix.txt",
               output_dir=config['output_dir'], sample=config["sampleList"],count_soft=config["count_soft"]),
        expand("{output_dir}/final_log/{sample}.log",
               output_dir=config['output_dir'], sample=config["sampleList"]),
     

