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

from demo.blueprint.modules.utils_biology import MappingConfigForm, get_biology_custom
from demo.blueprint.modules.utils_common import generate_steps_config,get_samples_list

biology= Blueprint('biology', __name__)


@biology.route('/', methods=['GET', 'POST'])
def index():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('biology.home', task_id=task_id))

@biology.route('/home/<task_id>', methods=['GET', 'POST'])
def home(task_id):
    form = MappingConfigForm()
    if form.validate_on_submit() and 'submit' in request.form:
        return redirect(url_for('biology.check',  task_id=task_id), code=307)
    return render_template('biology/home.html', form=form, task_id=task_id)
#########################################################################################
@biology.route('/biology/get_preset', methods=['POST'])
def get_preset():
    data = request.json
    soft_list, steps = get_biology_custom(data.get("presetKey"))
    return jsonify({"soft_list":soft_list, "steps":steps})

@biology.route('/biology/check/<task_id>', methods=['GET', 'POST'])
def check(task_id):
    form = request.form.to_dict()
    session["task_id"]=task_id
    session["select_raw_folder"]=form["input_dir"]
    fastqs = [f"{f}_1.fq.gz {f}_2.fq.gz" if form["mode"] == "pair" else f"{f}.fq.gz" for f in get_samples_list(form["input_dir"], form["mode"])]
    df_data = pd.DataFrame(fastqs, index=range(1,len(fastqs)+1), columns=[f'total: {len(fastqs)}']).to_html(justify='center')
    return render_template('biology/check.html',form = form , df_data=df_data, task_id=task_id )

@biology.route('/biology/submit_config_script', methods=['POST'])
def submit_dir_script():
    data = request.json 
    select_raw_folder = session.get('select_raw_folder')
    task_id = session.get('task_id')
    soft_step = data.get('steps', [])
    samples_list = get_samples_list(select_raw_folder, soft_step[0]["mode"], "fq.gz")
    soft_config = current_app.config['soft_config']
    paths_dict = current_app.config['paths_dict']
    code_script = current_app.config['code_script']
    output_path =os.path.join(current_app.config['paths_dict']["flask_out"], task_id)
    os.makedirs(output_path, exist_ok=True)
    config_string, snakemake_string, my_code = generate_steps_config(soft_step, select_raw_folder, samples_list, task_id, soft_config, paths_dict, code_script )
    with open(os.path.join(output_path,"config.yaml"), "w") as file1, open( os.path.join(output_path,"snakefile"), "w") as file2:
        print(config_string, file=file1)
        print(snakemake_string, file=file2)
    return jsonify({'my_code': my_code, "config_string": config_string, "snakemake_rule": snakemake_string})  



