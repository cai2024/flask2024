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
