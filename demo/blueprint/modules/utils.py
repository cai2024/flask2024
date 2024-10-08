import os
import csv
import json
import pandas as pd
import numpy as np
from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm




def get_fq_list(indir, paired_end=True):
    fastqs = []
    if paired_end:
        for file in os.listdir(indir):
            if file.endswith('1.fq.gz') and file[:-7] + '2.fq.gz' in os.listdir(indir):
                fastqs.append(file[:-7] + "1.fq.gz" + " " + file[:-7] + "2.fq.gz")
            if file.endswith('1.fastq.gz') and file[:-10] + '2.fastq.gz' in os.listdir(indir):
                fastqs.append(file[:-10] + "1.fastq.gz" + " " + file[:-10] + "2.fastq.gz")
    else:
        fastqs = [f for f in os.listdir(indir) if f.endswith('.fq.gz') or f.endswith('.fastq.gz')]
    return fastqs


def get_samples_list(select_raw_folder, fq_type, raw_myOptions):
    samples_list = []
    for filename in os.listdir(select_raw_folder):
        base_name = os.path.basename(filename)
        if filename.endswith(".fq.gz") and fq_type == "single" and raw_myOptions == "raw.fq.gz":
            samples_list.append(base_name[:-6])
        elif filename.endswith("_1.fq.gz") and fq_type == "pair" and raw_myOptions == "raw.fq.gz" and (
                filename[:-8] + "_2.fq.gz" in os.listdir(select_raw_folder)):
            samples_list.append(base_name[:-8])
        elif filename.endswith(".bam") and raw_myOptions == "raw.bam":
            samples_list.append(base_name[:-4])
    return samples_list

# name 数据量G cleanread比例 去重比例  比对率 基因组cover 平均depth  gc mit q20 q30 totalreads



def ana_get_res(outdir):
    target_keys = ['name', 'sequence_size', 'clean_read_ratio', 'duplication_rate',
                   'map_read_ratio', 'cover_ratio', 'average_depth', 'gc_filter_ratio',
                   'mit_ratio', 'q20_filter_ratio', 'q30_filter_ratio', 'total_reads_raw']

    if not os.path.exists(outdir):
        empty_data = pd.DataFrame(columns=target_keys)
        return empty_data


    data_list = []
    res_list = [os.path.join(outdir, f) for f in os.listdir(outdir) if f.endswith('.log')]
    for res_file in res_list:
        with open(res_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            keys = next(reader)
            values = next(reader)
            values = [f'0{value}' if value.startswith('.') else value for value in values]
            my_dir = dict(zip(keys, values))
            data_list.append(my_dir)

    # 一次性从列表创建DataFrame
    data = pd.DataFrame(data_list)
    target_keys=keys
    # 重新排序列以确保它们与target_keys的顺序匹配
    data = data[target_keys]
    data['index_name'] = data['name']
    data = data.set_index('index_name')
    data.index = range(1, data.shape[0] + 1)
    

    csv_path = os.path.join(outdir, 'analysis_results.csv')
    data.to_csv(csv_path)

    return data

