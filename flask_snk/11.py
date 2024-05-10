from os.path import join
from os.path import exists


# 定义主要软件依赖
trim_galore_path = "/data/biosoft/soft2024/trim_galore"
bwa_path = "/data/biosoft/soft2024/bwa"

# 从config.yaml加载配置
# 由命令行定义


# 规则
rule all:
    input:
        expand("{trim_out_path}/{sample}_1_val_1.fq.gz",  
            sample=config["sampleList"],trim_out_path=config['trim_out_path']),



# 对原始数据进行trimming
rule trim_galore_pe: 
    input:
        expand("{fq_in_path}/{{sample}}_{readDirection}.fq.gz",
            fq_in_path=config["fq_in_path"],
            readDirection=['1','2'])
    output:
        expand("{trim_out_path}/{{sample}}_{readDirection}_val_{readDirection}.fq.gz",
            trim_out_path=config["trim_out_path"],
            readDirection=['1','2'])
    log:
        "logs/trim_galore/{sample}.log"
    params:
        out_path=config['trim_out_path'],
    conda:
        'flask2024',
    threads: 20
    shell:
        """
        (trim_galore --cores {threads} --gzip -o {params.out_path} --paired {input}) 2> {log}
        """
