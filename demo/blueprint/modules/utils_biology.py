import os
import csv
import json
import pandas as pd
import numpy as np
from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
import yaml
from collections import defaultdict
from .utils_common import SoftStep, SoftStepManager, process_all_options


biology_preset_list=["WGS","WGBS","CAPS","CS"]


def get_biology_simple(preset_type,mode):
    biology_preset_custom = {
        "WGS":[{"step0 trim" : ["Trim trim_galore", "Trim fastp", "Trim cutadapt"]},
               {"step1 map" : ["map bwa", "map bowtie2"]},
        ],
        "WGBS":[{"step0 trim" : ["Trim trim_galore", "Trim fastp", "Trim cutadapt"]},
               {"step1 map" : ["map bwa", "map bowtie2"]},
        ],
    }
    return biology_preset_custom[preset_type]  

def get_biology_preset(preset_type,mode,soft_config,index_options_dict):
    preset_dict = {
        "WGS": {
            "baseOptions": {"pair": {"fq1": ["raw_1.fq.gz"], "fq2": ["raw_2.fq.gz"]}, "single": {"fq": ["raw.fq.gz"]}},
            "pair": [
                SoftStep(0, "Trim fastp", { }, {}, "", "pair", "raw.fq.gz",{}, 
                         {"fq1": "raw_1.fq.gz", "fq2": "raw_2.fq.gz"}, {}),
                SoftStep(1, "map bwa", {}, {}, "", "pair", "raw.fq.gz", {},
                         {"fq1": "step0_fastp_raw_1.fq.gz", "fq2": "step0_trim_galore_raw_2.fq.gz"},{}), ],
            "single": [
                SoftStep(0, "Trim fastp", {}, {}, "", "single", "raw.fq.gz", {}, 
                         {"fq": "raw.fq.gz"}, {}),
                SoftStep(1, "map bwa", {}, {}, "", "single", "raw.fq.gz", {},
                         {"fq": "step0_fastp_raw.fq.gz"}, {}), ]
        },
    }
    return  init_div_biology( [step.to_json() for step in preset_dict.get(preset_type, preset_dict["WGS"])[mode]],soft_config,index_options_dict) 

def init_div_biology(steps,soft_config,index_options_dict):
    manager = SoftStepManager()
    baseOptions = defaultdict(list)

    if steps[0]["mode"] == "pair" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq1"] = ["raw_1.fq.gz"]
        baseOptions["fq2"] = ["raw_2.fq.gz"]
    elif steps[0]["mode"] == "single" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq"] = ["raw.fq.gz"]
    else:
        baseOptions["bam"] = ["raw.bam"]

    all_options_list = process_all_options(soft_config, steps, baseOptions)

    for i, step in enumerate(steps):
        all_options = all_options_list[i].copy()
        soft_step = SoftStep(**step)
        key_list = soft_config.get(step["type"].split()[1], soft_config[step["type"].split()[0]])[step["mode"]]["input"].keys()

        soft_step.update_details(
            selectableOptions={key: all_options.get(key.split("_")[0], []).copy() for key in key_list})
       
        if not soft_step.rec_option:
            soft_step.update_details(rec_option={
                key: all_options[key][
                    max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                               for temp in all_options[key]]), key=lambda x: x[1])[0]] for key in key_list})
        soft_step.update_details(
            selectableIndex=index_options_dict.get(step["type"].split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()

class MappingConfigForm(FlaskForm):
    output_dir=HiddenField()
    input_dir = StringField('输入文件夹', validators=[DataRequired()], description='请输入文件所在文件夹')
    mode = SelectField("选择测序数据模式", choices=[('pair', '双端测序'), ('single', '单端测序')],validate_choice=False, validators=[DataRequired()], description='')
    pipeline = SelectField("选择分析流程" ,choices = biology_preset_list ,validate_choice=False, validators=[DataRequired()], description='')
    param_choice = SelectField('选择参数类型',choices=[('default', 'default'), ('custom', 'custom')], default='default', validators=[DataRequired()])
    cpu_count = IntegerField('CPU数量', validators=[DataRequired()],default=20 ,description='CPU数量，必须为正整数，不超过机器最大CPU数量')
    submit = SubmitField('开始比对')




