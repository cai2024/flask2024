from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime

from demo.blueprint.modules.utils_common import  load_config,generate_steps_config,get_samples_list

individual_code = Blueprint('individual_code', __name__)
available_steps={}




@individual_code.route('/individual_code/', methods=['GET', 'POST'])
def individual_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('individual_code.code_for_batch', task_id=task_id))
@individual_code.route('/individual_code/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('individual_code/home.html', available_steps=available_steps, task_id=task_id)


@individual_code.route('/individual_code/get_preset', methods=['POST'])
def get_preset():
    data = request.json 
    return jsonify(batch_preset(data.get("presetType"),data.get("mode")))

@individual_code.route('/individual_code/init_div_steps', methods=['POST'])
def init_div_steps():
    data = request.json
    return jsonify(init_div_steps3(data.get("steps")))

