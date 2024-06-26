FROM continuumio/anaconda3

# 设置工作目录

WORKDIR /app

# 添加channels并设置优先级
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge



# 安装Snakemake
RUN conda install -c bioconda snakemake=8.11.1 -y


# 安装其他软件包
RUN conda create --name flask2024 python=3.9 -y && \
    /bin/bash -c "source activate flask2024 && \
    conda install -c bioconda seqtk -y && \
    conda install -c bioconda cutadapt -y && \
    conda install -c bioconda trim-galore -y && \
    conda install samtools -y && \
    pip install gunicorn click flask flask-wtf bootstrap-flask pandas numpy PyYAML flask_sqlalchemy Bio Levenshtein pysam"


COPY . /app


# 确保虚拟环境被激活
ENV PATH /opt/conda/envs/flask2024/bin:$PATH
RUN sed -i '/"flask_out"/c\    "flask_out": "/data/cailab/flask_out",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"ref_path"/c\    "ref_path": "/data/reference2024",' /app/demo/blueprint/modules/paths.py


RUN sed -i '/"snakemake"/c\    "snakemake": "/data/biosoft/soft2024/conda/anaconda_23.4.7/bin/snakemake222",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"snakemake"/c\    "snakemake": "/opt/conda/bin/snakemake",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"java"/c\    "java": "/opt/conda/envs/flask2024/bin/java",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"fastp"/c\    "fastp":  "/app/soft/fastp-0.23.4/fastp",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bowtie2"/c\    "bowtie2": "/app/soft/bowtie2-2.5.2-linux-x86_64/bowtie2",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bowtie2_path"/c\    "bowtie2_path": "/app/soft/bowtie2-2.5.2-linux-x86_64",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bwa"/c\    "bwa": "/app/soft/bwa-0.7.17/bwa",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bismark"/c\    "bismark": "/app/soft/Bismark-0.24.2/bismark",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bismark_methylation_extractor"/c\    "bismark_methylation_extractor": "/app/soft/Bismark-0.24.2/bismark_methylation_extractor",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"picard"/c\    "picard": "/app/soft/picard/picard.jar",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"qualimap"/c\    "qualimap": "/app/soft/qualimap_v2.3/qualimap",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bedtools"/c\    "bedtools": "/app/soft/bedtools2/bin/bedtools",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"py_ref"/c\    "py_ref": "/app/flask_snk/ref",' /app/demo/blueprint/modules/paths.py
RUN sed -i '/"bamUtil"/c\    "bamUtil": "/app/soft/bamUtil/bin/bam",' /app/demo/blueprint/modules/paths.py






RUN conda create -n py39 python=3.9 -y 
RUN /bin/bash -c "source activate py39 && \
    pip install astair"



# 启动容器时的默认命令
CMD ["tail", "-f", "/dev/null"]










