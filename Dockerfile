FROM continuumio/anaconda3

# 设置工作目录
WORKDIR /data/cailab/flask2024

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

# 确保虚拟环境被激活
ENV PATH /opt/conda/envs/flask2024/bin:$PATH

# 启动容器时的默认命令
CMD ["tail", "-f", "/dev/null"]
