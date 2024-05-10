import os
import multiprocessing
import json
import time
import pandas as pd
import numpy as np
import yaml
import subprocess
import uuid
from datetime import datetime
from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify
from flask_wtf import FlaskForm
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from wtforms.validators import DataRequired, Length, Optional, URL



from demo.blueprint.modules.utils import filter_folders
from demo.blueprint.modules.utils import get_fq_list
from demo.blueprint.modules.utils import ana_get_res
from demo.blueprint.modules.utils import MappingConfigForm
from demo.blueprint.modules.utils import MappingCheckForm
from demo.blueprint.modules.utils import ResultForm
from demo.blueprint.modules.snk_config import submit_snakemake
from demo.blueprint.modules.core import PipelineForm

pipelineForm=PipelineForm()
biology= Blueprint('biology', __name__)
user_dir="/data/cailab/flask_output/"



@biology.route('/', methods=['GET', 'POST'])
def index():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('biology.home', task_id=task_id))

@biology.route('/home/<task_id>', methods=['GET', 'POST'])
def home(task_id):
    form = MappingConfigForm()
    form.pipeline.choices=pipelineForm.suffixes
    if form.validate_on_submit() and form.submit.data:
        return redirect(url_for('biology.check',  task_id=task_id), code=307)
    return render_template('biology/home.html', form=form, task_id=task_id)

@biology.route('/get_folders')
def get_folders():
    return jsonify(filter_folders(request.args.get('q', '')))
@biology.route('/process_folder_selection', methods=['POST'])
def process_folder_selection():
    task_id = request.form.get('task_id')
    session[f'selected_folder_{task_id}'] = request.form['selected_folder']
    return jsonify({'message': 'Folder selection received successfully'})

@biology.route('/check/<task_id>', methods=['GET','POST'])
def check(task_id):
    info = json.loads(json.dumps(request.form))
    selected_folder = session.get(f'selected_folder_{task_id}', 'no input')  
    fastqs = get_fq_list(selected_folder, info["mode"] != "single")
    df_data = pd.DataFrame(fastqs, index=range(1,len(fastqs)+1), columns=[f'total: {len(fastqs)}']).to_html(justify='center')
    form = pipelineForm.create_instance(info['pipeline']+"_"+info['param_choice'])
    form.data_ana.data = selected_folder
    form.output_dir.data = user_dir + task_id
    form.cpu_count.data = info['cpu_count']
    form.mode.data=info["mode"]
    if form.submit_ana.data and form.validate_on_submit():
        if not os.path.exists(info['output_dir']):
            os.makedirs(info['output_dir'])
        back_process = multiprocessing.Process(name=os.path.split(info['output_dir'])[1], target=submit_snakemake, args=(info,))
        back_process.daemon = True
        back_process.start()
        return redirect(url_for('biology.query',  task_id=task_id), code=307)
    return render_template('biology/check.html', data_info=df_data, task_id=task_id,  form=form)






@biology.route('/query/<task_id>', methods=['POST' , 'GET'])
def query(task_id):
    upload_form=ResultForm()
    return render_template('biology/query.html',upload_form=upload_form, task_id=task_id)

@biology.route('/fetch-logs')
def fetch_logs():
    task_id = request.args.get('task_id')
    with open( user_dir+ task_id +"/snk.log", 'r') as file:
        return jsonify({'logs': file.read()})

@biology.route('/progress/', methods=['POST'])
def get_progress():
    task_id = request.form.get('task_id')
    outdir=user_dir +task_id  + "/final_log"
    data_info = ana_get_res(outdir)
    return [ data_info.to_html()]


