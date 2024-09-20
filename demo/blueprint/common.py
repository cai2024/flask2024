from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context, g, \
    session, Flask, jsonify, send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from .modules.paths import paths_dict
from demo.blueprint.modules.utils import get_samples_list
from demo.blueprint.modules.utils_complex import generate_steps_config, generate_present_script
from demo.blueprint.modules.utils_complex import available_steps, my_code_simple, init_div_steps2






common_blueprint = Blueprint('common_blueprint', __name__)

@common_blueprint.route('/common/init_type_div', methods=['POST'])
def init_type_div():
    data = request.json
    result = get_init_div(data.get('step_id'), data.get('step_type'), data.get('step_params'), data.get('options'),
                          data.get('mode'))
    return jsonify({'div_string': result[0], 'rec_option': result[1]})

@common_blueprint.route('/common/update_my_options', methods=['POST'])
def update_my_options():
    data = request.json
    temp_options = get_result_option(str(data.get('step_id')), data.get('step_type'), data.get('mode'),
                                     data.get('temp_options'))
    return jsonify({'temp_options': {key: list(set(value)) for key, value in temp_options.items()}})

@common_blueprint.route('/common/index_select', methods=['POST'])
def index_select():
    data = request.json
    return jsonify({'index_options': get_type_index(data.get('stepId'), data.get('step_type'))})

@common_blueprint.route('/common/update_pair_select', methods=['POST'])
def update_pair_select():
    data = request.json
    return jsonify({'update_select_dir': update_select_dir(data.get('key_dir'), data.get('change_key'),
                                                           data.get('step_type'), data.get('mode'))})

@common_blueprint.route('/common/get_display_value', methods=['POST'])
def get_display_value():
    data = request.get_json()
    displayString = my_code_simple[data.get('step_type')][data.get('mode')].format(params=data.get('params'),
                                                                                   stepId=data.get('stepId'))
    return jsonify({'displayString': displayString})

@common_blueprint.route('/common/get_preset', methods=['POST'])
def get_preset():
    data = request.json
    presetType = data.get("presetType")
    mode = data.get("mode")
    if presetType == "WGS" and mode == "single":
        options = ['Trim trim_galore', 'Trim fastp', 'Trim trim_galore']
        select = [{'fq': 'raw.fq.gz'}, {'fq': 'step0_trim_galore.raw.fq.gz'}, {'fq': 'step1_fastp.raw.fq.gz'}]
        params = ["-A NNNN", "", "-D"]

    return jsonify({"options": options, "select": select, "params": params})

@common_blueprint.route('/common/submit_form_steps', methods=['POST'])
def submit_form_steps():
    data = request.json
    session["fq_type"] = data.get("fq_type")
    session["raw_myOptions"] = data.get("raw_myOptions")
    session["steps"] = data.get("steps")
    return jsonify({"step": data.get("steps")})





@common_blueprint.route('/common/init_div_steps', methods=['POST'])
def init_div_steps():
    data = request.json
    return jsonify(init_div_steps2(data.get("steps")))

