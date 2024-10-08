import os
import csv
import json
import pandas as pd
import numpy as np
from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm

biology_preset_list=["WGS","WGBS","CAPS","CS"]



class MappingConfigForm(FlaskForm):
    output_dir=HiddenField()
    input_dir = StringField('输入文件夹', validators=[DataRequired()], description='请输入文件所在文件夹')
    mode = SelectField("选择测序数据模式", choices=[('pair', '双端测序'), ('single', '单端测序')],validate_choice=False, validators=[DataRequired()], description='')
    pipeline = SelectField("选择分析流程" ,choices = biology_preset_list ,validate_choice=False, validators=[DataRequired()], description='')
    param_choice = SelectField('选择参数类型',choices=[('default', 'default'), ('custom', 'custom')], default='default', validators=[DataRequired()])
    cpu_count = IntegerField('CPU数量', validators=[DataRequired()],default=20 ,description='CPU数量，必须为正整数，不超过机器最大CPU数量')

    submit = SubmitField('开始比对')




