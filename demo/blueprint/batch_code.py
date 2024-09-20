from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from .modules.paths import paths_dict
from demo.blueprint.modules.utils import get_samples_list
from demo.blueprint.modules.utils_complex import generate_steps_config,generate_present_script, batch_preset
from demo.blueprint.modules.utils_complex import available_steps, my_code_simple

complex_code = Blueprint('complex_code', __name__)
# 分离出实际有效的公共路由

@complex_code.route('/complex_code/', methods=['GET', 'POST'])
def code_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('complex_code.code_for_batch', task_id=task_id))
@complex_code.route('/code_for_batch/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('complex_code/home.html', available_steps=available_steps, task_id=task_id)


@complex_code.route('/batch_code/get_preset', methods=['POST'])
def get_preset():
    data = request.json 
    return jsonify(batch_preset(data.get("presetType"),data.get("mode")))



################################################################################################################


@complex_code.route('/batch_code/store_soft_step', methods=['POST'])
def store_soft_step():
    data = request.json
    session['soft_step'] = data['soft_step']
    return '', 204

@complex_code.route('/batch_code/submit_form_code/<task_id>', methods=['GET', 'POST'])
def submit_form_code(task_id):
    session["task_id"] = task_id
    return render_template('complex_code/code.html', task_id=task_id)

#####################################



@complex_code.route('/batch_code/submit_config_script', methods=['POST'])
def submit_dir_script():

    data = request.json
    select_raw_folder = data.get('select_raw_folder')
    task_id = data.get("task_id")
    soft_step = session.get('soft_step', [])
    samples_list = get_samples_list(select_raw_folder, soft_step[0]["mode"], soft_step[0]["raw_type"])
    

    config_string, snakemake_string = generate_steps_config(soft_step, select_raw_folder, soft_step[0]["mode"], task_id, samples_list, paths_dict)
    my_code = generate_present_script(soft_step, soft_step[0]["mode"], select_raw_folder, samples_list, task_id, paths_dict)


    os.makedirs("/data/cailab/flask_out/" + task_id, exist_ok=True)
    with open("/data/cailab/flask_out/" + task_id + "/config.yaml", "w") as file1, open( "/data/cailab/flask_out/" + task_id + "/snakefile", "w") as file2:
        print(config_string, file=file1)
        print(snakemake_string, file=file2)
    return jsonify({'my_code': my_code, "config_string": config_string, "snakemake_rule": snakemake_string})
   


@complex_code.route('/submit_smk_job', methods=['POST'])
def submit_smk_job():
    data = request.json
    task_id=data.get('task_id')
    outdir="/data/cailab/flask_out/" + task_id
    snakemake_cmd = [
        'nohup',
        'snakemake',
        '--snakefile', outdir+'/snakefile' ,
        '-c', '10',
        '--use-conda',
        '--conda-frontend', 'conda',
        '--configfile', outdir+'/config.yaml',
        '--default-resources' ,  f"tmpdir='{outdir}/temp_test'", 'mem_mb=200000', 'disk_mb=200000', 
       
    ]

    try:
        # 使用Popen使命令在后台运行
        with open(outdir + '/snk.log', 'w') as f:
            subprocess.Popen(snakemake_cmd, stdout=f, stderr=subprocess.STDOUT, preexec_fn=os.setpgrp)
        # 立即返回任务ID，客户端可用task_id来查询任务状态
        return jsonify({'task_id': task_id })
    except Exception as e:
        with open(outdir + '/snk.log', 'a') as f:
            print(f"运行出错: {e}", file=f)
            print("请确认参数是否正确，或者@wangzc", file=f)
        return jsonify({'error': '任务启动失败，请检查日志', 'task_id': task_id}), 500


