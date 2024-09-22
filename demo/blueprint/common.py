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




@common_blueprint.route('/common/get_display_value', methods=['POST'])
def get_display_value():
    data = request.get_json()
    displayString = my_code_simple[data.get('step_type')][data.get('mode')].format(params=data.get('params'),
                                                                                   stepId=data.get('stepId'))
    return jsonify({'displayString': displayString})




