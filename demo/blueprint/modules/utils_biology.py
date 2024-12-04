import os
import csv
import json
import pandas as pd
import numpy as np
import yaml
import copy
import gzip
import base64
from io import BytesIO
from graphviz import Digraph
from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from collections import defaultdict
from .utils_common import SoftStep, SoftStepManager, process_all_options


biology_preset_list=["WGS","WGBS","CAPS","CS"]
class MappingConfigForm(FlaskForm):
    output_dir=HiddenField()
    input_dir = StringField('输入文件夹', validators=[DataRequired()], description='请输入文件所在文件夹')
    mode = SelectField("选择测序数据模式", choices=[('pair', '双端测序'), ('single', '单端测序')],validate_choice=False, validators=[DataRequired()], description='')
    pipeline = SelectField("选择分析流程" ,choices = biology_preset_list ,validate_choice=False, validators=[DataRequired()], description='')
    param_choice = SelectField('选择参数类型',choices=[('default', 'default'), ('custom', 'custom')], default='default', validators=[DataRequired()])
    cpu_count = IntegerField('CPU数量', validators=[DataRequired()],default=20 ,description='CPU数量，必须为正整数，不超过机器最大CPU数量')
    submit = SubmitField('开始比对')

def get_biology_custom(base64_encoded_data):
    categories = [
        ["Trim trim_galore", "Trim fastp", "Trim cutadapt"],
        ["Trim hairpin_cut", "Trim hairpin_cut_N"],
        ["map bwa", "map bowtie2", "map bismark", "map star"],
        ["map bwa_spikein", "map bowtie2_spikein", "map bismark_spikein"],
    ]
    category_map = {item: group for group in categories for item in group}
    compressed_data = base64.urlsafe_b64decode(base64_encoded_data.encode('utf-8'))
    json_data = gzip.decompress(compressed_data).decode('utf-8')
    original_data = json.loads(json_data)
    soft_list=[]
    manager = SoftStepManager()
    for i,temp in enumerate(original_data):
        soft_list.append({f"step{i} {temp[1].split()[0]}": category_map.get(temp[1],[temp[1]])})  
        manager.add_step(SoftStep(temp[0],temp[1],{},{},temp[2],temp[3],temp[4],{},temp[5],{})) 
    return  soft_list, manager.to_json()

