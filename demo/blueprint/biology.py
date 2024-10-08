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



from demo.blueprint.modules.utils_biology import MappingConfigForm
from demo.blueprint.modules.utils_complex import generate_steps_config,generate_present_script, batch_preset, init_div_steps2
from demo.blueprint.modules.utils_complex import available_steps, my_code_simple

from demo.blueprint.modules.utils import get_fq_list
from demo.blueprint.modules.utils import ana_get_res
from demo.blueprint.modules.snk_config import submit_snakemake
from demo.blueprint.modules.paths import paths_dict

biology= Blueprint('biology', __name__)
user_dir=paths_dict["flask_out"]



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

@biology.route('/biology/check/<task_id>', methods=['GET', 'POST'])
def check(task_id):
   
    form = request.form.to_dict()
    input_dir=form["input_dir"]
    mode=form["mode"]
    pipeline=form["pipeline"]

    fastqs = get_fq_list(form["input_dir"], form["mode"] != "single")
    df_data = pd.DataFrame(fastqs, index=range(1,len(fastqs)+1), columns=[f'total: {len(fastqs)}']).to_html(justify='center')
    return render_template('biology/check.html',form = form , df_data=df_data, task_id=task_id )


biology_preset_dict = {
    "WGS":[{"step":"step0 trim","select":["Trim trim_galore", "Trim_fastp", "Trim_cutadapt"]},
           {"step":"step1 map","select":["map bwa", "map bowtie2"]},
    ],
    "WGBS":[{"step":"step0 trim","select":["Trim trim_galore", "Trim_fastp", "Trim_cutadapt"]},
           {"step":"step1 map","select":["map bwa", "map bowtie2"]},
    ],    
}

@biology.route('/biology/get_preset', methods=['POST'])
def get_preset():
    data = request.json 
    return jsonify(biology_preset_dict[data.get("form")["pipeline"]])



#########################################################################################
@biology.route('/biology/submit_form_code/<task_id>', methods=['GET', 'POST'])
def submit_form_code(task_id):
    if request.method == 'POST':
        data = request.json
        session["task_id"] = task_id
        session['soft_step'] = data['soft_step']
        return jsonify({"task_id": task_id}) 
    elif request.method == 'GET':
        return render_template('complex_code/code.html', task_id=task_id)



@biology.route('/biology/submit_config_script', methods=['POST'])
def submit_dir_script():
    data = request.json
    select_raw_folder = data.get('select_raw_folder')
    task_id = data.get("task_id")
    soft_step = session.get('soft_step', [])
    samples_list = get_samples_list(select_raw_folder, soft_step[0]["mode"], soft_step[0]["raw_type"])

    config_string, snakemake_string = generate_steps_config(soft_step, select_raw_folder , samples_list, task_id, paths_dict)
    my_code = generate_present_script(soft_step, select_raw_folder, samples_list, task_id, paths_dict)
    with open("/data/cailab/flask_out/" + task_id + "/config.yaml", "w") as file1, open( "/data/cailab/flask_out/" + task_id + "/snakefile", "w") as file2:
        print(config_string, file=file1)
        print(snakemake_string, file=file2)
    return jsonify({'my_code': my_code, "config_string": config_string, "snakemake_rule": snakemake_string})


#########################################################################################
@biology.route('/query/<task_id>', methods=['POST' , 'GET'])
def query(task_id):
    upload_form=MappingConfigForm()
    return render_template('biology/query.html',upload_form=upload_form, task_id=task_id)

@biology.route('/fetch-logs')
def fetch_logs():
    task_id = request.args.get('task_id')
    with open( user_dir+"/" +task_id +"/snk.log", 'r') as file:
        return jsonify({'logs': file.read()})

@biology.route('/progress/', methods=['POST'])
def get_progress():
    task_id = request.form.get('task_id')
    outdir=user_dir+"/" +task_id  + "/final_log"
    data_info = ana_get_res(outdir)
    return [ data_info.to_html()]


