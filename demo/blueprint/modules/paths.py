# 注意路径末尾不要带/

paths_dict = {
    # soft路径
    "snakemake": "/data/biosoft/soft2024/conda/anaconda_23.4.7/bin/snakemake", #容器conda额外修改 
    # trim软件
    "fastp":  "/data/cailab/flask2024/soft/fastp-0.23.4/fastp",
    "trim_galore": "trim_galore", #conda
    "cutadapt": "cutadapt", #conda
    # 比对软件
    "bowtie2": "/data/cailab/flask2024/soft/bowtie2-2.5.2-linux-x86_64/bowtie2",
    "bowtie2_path": "/data/cailab/flask2024/soft/bowtie2-2.5.2-linux-x86_64",
    "bwa": "/data/cailab/flask2024/soft/bwa-0.7.17/bwa",
    "bismark": "/data/cailab/flask2024/soft/Bismark-0.24.2/bismark",
    "bismark_methylation_extractor": "/data/cailab/flask2024/soft/Bismark-0.24.2/bismark_methylation_extractor",
    "STAR" : "/data/cailab/flask2024/soft/STAR-2.7.11b/bin/Linux_x86_64_static/STAR",
    # other
    "picard": "/data/cailab/flask2024/soft/picard/picard.jar",
    "qualimap": "/data/cailab/flask2024/soft/qualimap_v2.3/qualimap",
    "bamUtil": "/data/cailab/flask2024/soft/bamUtil/bin/bam",
    "featureCounts" : "/data/cailab/flask2024/soft/subread-2.0.6-Linux-x86_64/bin/featureCounts",
    "bedtools": "/data/cailab/flask2024/soft/bedtools2/bin/bedtools",
    "samtools": "samtools",
    "java": "/usr/bin/java", #存疑
   # 配置目录
    "flask_out": "/data/cailab/flask_out" ,
    "ref_path": "/data/reference2024",
    "py_ref": "/data/cailab/flask2024/flask_snk/ref",   
}


"""
paths_dict = {
    "snakemake": "/home/soft/conda/anaconda_23.4.7/bin/snakemake",
    "fastp": "/home/soft/fastp/fastp-0.23.4/fastp",
    "bowtie2": "/home/soft/bowtie2/bowtie2-2.5.2-linux-x86_64/bowtie2",
    "bowtie2_path": "/home/soft/bowtie2/bowtie2-2.5.2-linux-x86_64",
    "bwa": "/home/soft/bwa/bwa-0.7.17/bwa",
    "bismark": "/home/soft/bismark/Bismark-0.24.2/bismark",
    "bismark_methylation_extractor": "/home/soft/bismark/Bismark-0.24.2/bismark_methylation_extractor",
    "STAR" : "/home/soft/star/STAR-2.7.11b/bin/Linux_x86_64_static/STAR",
    "picard": "/home/soft/picard/picard.jar",
    "qualimap": "/home/soft/qualimap/qualimap_v2.3/qualimap",
    "bamUtil": "/home/soft/bamUtil/bin/bam",
    "featureCounts": "/home/soft/subread/subread-2.0.6-Linux-x86_64/bin",
    "bedtools": "/home/soft/bedtools/bedtools2/bin/bedtools",
    "java": "/usr/bin/java",
    "ref_path": "/home/reference",
    "py_ref": "/home/flask/flask2024/flask_snk/ref",
    "flask_out": "/home/flask/flask_out",
}
"""
