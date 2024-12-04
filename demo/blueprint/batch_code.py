from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort, g
import os
import redis
import subprocess
import yaml
import uuid
import json
from datetime import datetime
from demo.blueprint.modules.utils_batch_code import get_batch_preset, init_div_batch, look_div_batch
from demo.blueprint.modules.utils_common import generate_steps_config, get_samples_list
batch_code = Blueprint('batch_code', __name__)
# 分离出实际有效的公共路由




@batch_code.route('/batch_code/', methods=['GET', 'POST'])
def batch_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('batch_code.code_for_batch', task_id=task_id))
@batch_code.route('/code_for_batch/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('batch_code/home.html', available_steps=current_app.config['available_steps'],available_steps2=current_app.config['available_steps2'] , task_id=task_id)

@batch_code.route('/batch_code/get_preset', methods=['POST'])
def get_preset():
    data = request.json
    return jsonify(get_batch_preset(data.get("base64_encoded_data"),current_app.config['soft_config'],current_app.config['index_options_dict'])) 

@batch_code.route('/batch_code/init_div_steps', methods=['POST'])
def init_div_steps():
    data = request.json
    return jsonify(init_div_batch(data.get("steps"),current_app.config['soft_config'],current_app.config['index_options_dict']))

@batch_code.route('/batch_code/display_div_steps', methods=['POST'])
def display_div_steps():
    data = request.json
    base64_encoded_data, img_base64=look_div_batch(data.get("steps"),current_app.config['soft_config'],current_app.config['index_options_dict'])
    return jsonify({"base64_encoded_data": base64_encoded_data, 'image': f"data:image/png;base64,{img_base64}"})
################################################################################################################

@batch_code.route('/batch_code/submit_form_code/<task_id>', methods=['GET', 'POST'])
def submit_form_code(task_id):
    if request.method == 'POST':
        data = request.get_json()
        session['steps'] = data.get('steps')
        session["task_id"] = task_id
        return jsonify({"task_id": task_id}) 
    elif request.method == 'GET':
        return render_template('batch_code/present_code.html', task_id=task_id)

@batch_code.route('/batch_code/submit_config_script', methods=['POST'])
def submit_dir_script():
    data = request.get_json()
    select_raw_folder = data.get('select_raw_folder')
    task_id = session.get("task_id")
    steps = session.get('steps', [])
    samples_list = get_samples_list(select_raw_folder, steps[0]["mode"] , list(steps[0]["baseOptions"].values())[0][0][4:]  )
    soft_config = current_app.config['soft_config']
    paths_dict = current_app.config['paths_dict']
    code_script = current_app.config['code_script']

    output_path =os.path.join(current_app.config['paths_dict']["flask_out"], task_id)
    os.makedirs(output_path, exist_ok=True)

    config_string, snakemake_string, my_code = generate_steps_config(steps, select_raw_folder, samples_list, task_id, soft_config, paths_dict, code_script )
    output_path =os.path.join(current_app.config['paths_dict']["flask_out"], task_id)
    with open(os.path.join(output_path,"config.yaml"), "w") as file1, open( os.path.join(output_path,"snakefile"), "w") as file2:
        print(config_string, file=file1)
        print(snakemake_string, file=file2)
    return jsonify({'my_code': my_code, "config_string": config_string, "snakemake_rule": snakemake_string})   


