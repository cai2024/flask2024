from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from .modules.paths import paths_dict
from demo.blueprint.modules.utils_complex import get_result_option, get_type_option, get_samples_list, generate_steps_config,generate_present_script
from demo.blueprint.modules.utils_complex import available_steps, soft_format

complex_code = Blueprint('complex_code', __name__)


@complex_code.route('/complex_code/', methods=['GET', 'POST'])
def code_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('complex_code.code_for_batch', task_id=task_id))

@complex_code.route('/code_for_batch/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('complex_code/home.html', available_steps=available_steps, task_id=task_id)


#####################################
@complex_code.route('/update_my_options', methods=['POST'])
def update_my_options():
    data = request.json
    my_options = data.get('my_options')
    my_options=get_result_option(str(data.get('step_id')),data.get('step_type'),data.get('select_option'),my_options)  # 新增选项保存到总选项
    return jsonify({'my_options': list(set(my_options))})

@complex_code.route('/index_select', methods=['POST'])
def update_index():
    data = request.json
    step_type = data.get('step_type').split()[1]
    index_options=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"]
    return jsonify({'index_options': index_options})

@complex_code.route('/get_options', methods=['POST'])
def get_options():
    data = request.json
    step_type = data.get('step_type')
    my_options = data.get('my_options')
    options = get_type_option(step_type,my_options)  # 公用的获取选项框应该提供的选项，也就是总选项的子集
    return jsonify({'options': options})

####################################
@complex_code.route('/get_display_value', methods=['POST'])
def update_display_params():
    data = request.get_json()
    step_type=data['type'].split()[1]
    step_params=data['params']
    select_option = data['select_option']
    if select_option[-3:]==".fq":
        input_file=select_option[:-3]
    else:
        input_file=select_option[:-4]
    displayString = soft_format[step_type].format(params=step_params,input=input_file)
    return jsonify({'displayString': displayString})

@complex_code.route('/get_preset', methods=['POST'])
def preset_soft():
    data = request.json
    presetType=data.get("presetType")
    if presetType=="WGS":
        options=['Trim trim_galore','map bowtie2', 'BAM picard']
        select=['raw.fq','step0_raw.fq', 'step1_step0_raw.fq']
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
    return render_template('complex_code/code.html', task_id=task_id)

# 从这里进入submit_form_code，也就是code.html
#####################################



@complex_code.route('/submit_dir_script', methods=['POST'])
def submit_dir_script():
    steps = session['steps']
    fq_type = session["fq_type"]
    raw_myOptions = session["raw_myOptions"]
    data = request.json
    task_id = data.get('task_id')
    select_raw_folder = data.get('select_raw_folder')

    samples_list = get_samples_list(select_raw_folder, fq_type, raw_myOptions)

    input_list = [record["select_option"] for record in steps]
    soft_list = [record["type"].split()[1] for record in steps]
    type_list = [record["type"].split()[0] for record in steps]
    params_list= [record["params"] for record in steps]
    config_steps, output_config_list = generate_steps_config(soft_list, type_list, input_list, select_raw_folder,fq_type , raw_myOptions)


    # 首先获取step字符串
    config_string = (f"mode: {fq_type}\n"
                   f"output_dir: /data/cailab/flask_out/{task_id}\n")
    config_string += "samples:\n" + "\n".join(f"  - {sample}" for sample in samples_list) + "\n"
    config_string += "\n".join(f"step{i}_params : \"{params}\"" for i,params in enumerate(params_list)) + "\n"
    config_string += config_steps
    # 获取snakemake
    snakemake_string=(
            "configfile: \"config.yaml\"\n"
            "include: \"/data/cailab/flask2024/smk/trim.smk\"\n"
            "rule all:\n"
            "    input:\n")
    snakemake_string += "\n".join( f"        expand({temp}, sample=config[\"samples\"], output_dir=config[\"output_dir\"])," for temp in [record.split("\n")[0].split(":")[1].strip() for record in output_config_list]) + "\n"
 
    my_code = (f"mode=\"{fq_type}\"\n"
               f"output_dir=\"/data/output/{task_id}\"\n")
    my_code += '\n'.join([f'{soft}="{paths_dict[soft]}"' for soft in soft_list]) + "\n"
    my_code += 'samples=("{}")'.format('" "'.join(samples_list)) + "\n"
    my_code += 'mkdir -p $output_dir/logs\n'
    my_code += 'for sample in "${samples[@]}"; do\n'
    my_code += generate_present_script([record["type"] for record in steps], fq_type, params_list, input_list, select_raw_folder) + "\n"
    my_code += "done\n"

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
        with open(outdir + '/snk.log', 'w') as f:
            subprocess.run(snakemake_cmd , stdout=f, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        with open(outdir + '/snk.log', 'a') as f:
            print(f"运行出错:{e}",file=f)
            print("请确认参数是否正确，或者@wangzc", file = f)
    return jsonify({'task_id': task_id })

