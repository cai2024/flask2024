from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import subprocess
import yaml
import uuid
from datetime import datetime
from .modules.paths import paths_dict
from demo.blueprint.modules.utils import get_samples_list
from demo.blueprint.modules.utils_complex import generate_steps_config,generate_present_script
from demo.blueprint.modules.utils_complex import available_steps, my_code_simple

individual_code = Blueprint('individual_code', __name__)


@individual_code.route('/individual_code/', methods=['GET', 'POST'])
def individual_home():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    return redirect(url_for('individual_code.code_for_batch', task_id=task_id))
@individual_code.route('/individual_code/<task_id>', methods=['GET', 'POST'])
def code_for_batch(task_id):
    return render_template('individual_code/home.html', available_steps=available_steps, task_id=task_id)



