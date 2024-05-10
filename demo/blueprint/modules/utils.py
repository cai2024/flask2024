import os
import csv
import json
import pandas as pd
import numpy as np
from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm


class MappingConfigForm(FlaskForm):
    output_dir=HiddenField()
    mode = SelectField("选择测序数据模式", choices=[('pair', '双端测序'), ('single', '单端测序')],validate_choice=False, validators=[DataRequired()], description='')
    pipeline = SelectField("选择分析流程" ,validate_choice=False, validators=[DataRequired()], description='')
    param_choice = SelectField('选择参数类型',
                              choices=[('default', 'default'), ('custom', 'custom')],
                              default='default', validators=[DataRequired()])
    cpu_count = IntegerField('CPU数量', validators=[DataRequired()],default=20 ,description='CPU数量，必须为正整数，不超过机器最大CPU数量')
    

    submit = SubmitField('开始比对')


class MappingCheckForm(FlaskForm):
    data_ana = StringField('数据目录', render_kw={'readonly':True})
    output_dir = StringField('输出目录', render_kw={'readonly':True})
    mode = HiddenField()
    pipeline= HiddenField()
    param_choice=HiddenField()
    cpu_count = HiddenField()


class ResultForm(FlaskForm):
    res_dir = StringField('结果目录', render_kw={'readonly':True})
    upload_res = SubmitField('上传结果', id='upload', render_kw={'hidden':True})









def filter_folders(query):
    folder_path = "/" if query.rfind("/") == 0 else query[: query.rfind("/")]
    folders = []
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        try:
            temp = [os.path.join(folder_path, d) for d in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, d))]
            folders.extend(temp)
            for i in temp:
                try:
                    subfolders = [os.path.join(i, subitem) for subitem in os.listdir(i) if os.path.isdir(os.path.join(i, subitem))]
                    folders.extend(subfolders)
                except PermissionError:
                    continue
        except PermissionError:
            pass
    query_level = query.count('/')
    matched_folders = [folder for folder in folders if folder.count('/') == query_level and query.lower() in folder.lower()]
    if len(matched_folders) == 1:
            matched_folders = [folder for folder in folders if (folder.count('/')-1  == query_level or folder.count('/') == query_level) and query.lower() in folder.lower()]
            return matched_folders
    return matched_folders



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

