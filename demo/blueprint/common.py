from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context, g, \
    session, Flask, jsonify, send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from demo.blueprint.modules.utils_common import generate_steps_config, get_samples_list 
import pandas as pd

common_blueprint = Blueprint('common_blueprint', __name__)


@common_blueprint.route('/common/submit_smk_job', methods=['POST'])
def submit_smk_job():
    data = request.json
    task_id = data.get('task_id')
    outdir = os.path.join(current_app.config['paths_dict']["flask_out"], task_id)
    snakemake_cmd = [
        'nohup',
        'snakemake',
        '--snakefile', outdir + '/snakefile',
        '-c', '10',
        '--use-conda',
        '--conda-frontend', 'conda',
        '--configfile', outdir + '/config.yaml',
        '--default-resources', f"tmpdir='{outdir}/temp_test'", 'mem_mb=200000', 'disk_mb=200000',
    ]
    try:
        with open(outdir + '/snk.log', 'w') as f:
            subprocess.Popen(snakemake_cmd, stdout=f, stderr=subprocess.STDOUT, preexec_fn=os.setpgrp)
        return jsonify({'task_id': task_id})
    except Exception as e:
        with open(outdir + '/snk.log', 'a') as f:
            print(f"运行出错: {e}", file=f)
            print("请确认参数是否正确，或者@wangzc", file=f)
        return jsonify({'error': '任务启动失败，请检查日志', 'task_id': task_id}), 500






@common_blueprint.route('/query/<task_id>', methods=['POST' , 'GET'])
def query(task_id):
    upload_form={}
    return render_template('common/query.html',upload_form=upload_form, task_id=task_id)

@common_blueprint.route('/common/fetchLogs',methods=['POST'])
def fetch_logs():
    data = request.get_json()
    with open( current_app.config['paths_dict']["flask_out"]+"/" +data.get("task_id") +"/snk.log", 'r') as file:
        return jsonify({'logs': file.read()})

@common_blueprint.route('/common/progress/', methods=['POST'])
def get_progress():
    data = request.get_json()
    output_dir = os.path.join(current_app.config['paths_dict']["flask_out"], data.get("task_id"))
    df = get_bate_table(output_dir)
    return [ df.to_html()]

def get_bate_table(output_dir,):
    config_data = yaml.safe_load(open(os.path.join(output_dir, "config.yaml")))
    result_dict = {step: {} for step in config_data["steps"]}
    sample_list = config_data["samples"]
    log_paths = [os.path.join(output_dir, f"logs/{sample}.log") for sample in sample_list if
                 os.path.exists(os.path.join(output_dir, f"logs/{sample}.log"))]

    for log_path in log_paths:
        with open(log_path) as file:
            for line in file:
                line=line.strip().split()
                if line[1] not in result_dict[line[0]]:
                    result_dict[line[0]][line[1]] = {}
                result_dict[line[0]][line[1]][os.path.splitext(os.path.basename(log_path))[0]] = line[2]
    table = []
    for step, step_dict in result_dict.items():
        for key, value in step_dict.items():
            table.append([step, key] + [value.get(sample) for sample in sorted(sample_list)])

    df = pd.DataFrame(table, columns=['Step', 'Param'] + sorted(sample_list)).T.reset_index()
    df.columns = df.iloc[0]
    df = df.drop(0)
    df.to_csv(os.path.join(output_dir, 'logs/analysis_results.csv'))
    return df






@common_blueprint.route('/common/get_display_value', methods=['POST'])
def get_display_value():
    data = request.get_json()
    displayString = my_code_simple[data.get('step_type')][data.get('mode')].format(params=data.get('params'),
                                                                                   stepId=data.get('stepId'))
    return jsonify({'displayString': displayString})



@common_blueprint.route('/common/get_folders')
def get_folders():
    return jsonify(filter_folders(request.args.get('folder', '')))

@common_blueprint.route('/common/store_folder', methods=['POST'])
def store_folder():
    data = request.get_json()
    session['selected_folder'] = data.get('folder')
    return jsonify({'message': 'Folder selection received successfully'})



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
