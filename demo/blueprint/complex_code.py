from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from .modules.paths import paths_dict
from demo.blueprint.modules.utils import get_samples_list
from demo.blueprint.modules.utils_complex import get_result_option, generate_steps_config,generate_present_script,get_type_index,get_init_div,update_select_dir
from demo.blueprint.modules.utils_complex import available_steps, my_code_simple

complex_code = Blueprint('complex_code', __name__)


@complex_code.route('/complex_code/', methods=['GET', 'POST'])
def code_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('complex_code.code_for_batch', task_id=task_id))
@complex_code.route('/code_for_batch/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('complex_code/home.html', available_steps=available_steps, task_id=task_id)


@complex_code.route('/init_type_div', methods=['POST'])
def init_type_div():
    data = request.json
    result= get_init_div(data.get('step_id'),data.get('step_type'),data.get('step_params'),data.get('options'),data.get('mode'))
    return jsonify({'div_string':result[0],'rec_option':result[1] }) 
@complex_code.route('/update_my_options', methods=['POST'])
def update_my_options():
    data = request.json
    temp_options=get_result_option(str(data.get('step_id')),data.get('step_type'),data.get('mode'),data.get('temp_options'))  
    return jsonify({'temp_options': {key: list(set(value)) for key, value in temp_options.items()} })
@complex_code.route('/index_select', methods=['POST'])
def update_index():
    data = request.json
    return jsonify({'index_options': get_type_index(data.get('stepId'),data.get('step_type'))})

################################################################################################################
@complex_code.route('/update_pair_select', methods=['POST'])
def update_pair_select():
    data = request.json
    return jsonify({'update_select_dir': update_select_dir(data.get('key_dir'),data.get('change_key'), data.get('step_type'),data.get('mode')) })

@complex_code.route('/get_display_value', methods=['POST'])
def get_display_value():
    data = request.get_json()
    displayString = my_code_simple[data.get('step_type')][data.get('mode')].format(params=data.get('params'),stepId=data.get('stepId'))
    return jsonify({'displayString': displayString})

##################################################################################################################

@complex_code.route('/get_preset', methods=['POST'])
def preset_soft():
    data = request.json
    presetType=data.get("presetType")
    mode=data.get("mode")
    if presetType=="WGS" and mode == "single":
        options=['Trim trim_galore','Trim fastp', 'Trim trim_galore']
        select=[{'fq':'raw.fq.gz'},{'fq':'step0_trim_galore.raw.fq.gz'}, {'fq':'step1_fastp.raw.fq.gz'}]
        params=["-A NNNN","","-D"]
    
    return jsonify({"options":options,"select":select,"params":params})



@complex_code.route('/submit_form_steps', methods=['POST'])
def submit_form_steps():
    data = request.json
    session["fq_type"] = data.get("fq_type")
    session["raw_myOptions"] = data.get("raw_myOptions")
    session["steps"] = data.get("steps")
    return jsonify({"step":data.get("steps")})

@complex_code.route('/submit_form_code/<task_id>', methods=['GET', 'POST'])
def submit_form_code(task_id):
    session["task_id"] = task_id
    return render_template('complex_code/code.html', task_id=task_id)

#####################################



@complex_code.route('/submit_dir_script', methods=['POST'])
def submit_dir_script():
    steps = session['steps']
    mode = session["fq_type"]
    raw_myOptions = session["raw_myOptions"]
    task_id = session["task_id"]
    data = request.json
    select_raw_folder = data.get('select_raw_folder')
    samples_list = get_samples_list(select_raw_folder, mode, raw_myOptions)

    config_string, snakemake_string = generate_steps_config(steps, select_raw_folder, mode, task_id, samples_list,paths_dict)
    my_code = generate_present_script(steps, mode, select_raw_folder, samples_list, task_id, paths_dict)


    os.makedirs("/data/cailab/flask_out/" + task_id, exist_ok=True)
    with open("/data/cailab/flask_out/" + task_id + "/config.yaml", "w") as file1, open(
            "/data/cailab/flask_out/" + task_id + "/snakefile", "w") as file2:
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


