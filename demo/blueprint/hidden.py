import os
import multiprocessing
import json
import time
import pandas as pd
import numpy as np
import yaml
import subprocess
import uuid
from datetime import datetime
from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify
from flask_wtf import FlaskForm
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from wtforms.validators import DataRequired, Length, Optional, URL

hidden_bp = Blueprint('hidden', __name__, template_folder='templates')



@hidden_bp.route('/login', methods=['POST'])
def login():
    task_id = request.form['task_id']  # 从表单获取task_id
    password = request.form['password']
    if password == 'wangzc':
        session[f'logged_in_{task_id}'] = True
        return redirect(url_for('hidden.hidden', task_id=task_id))
    else:
        return "密码错误！", 401

@hidden_bp.route('/hidden/<task_id>')
def hidden(task_id):
    if not session.get(f'logged_in_{task_id}'):
        return redirect(url_for('biology.home', task_id=task_id))
    return render_template('hidden/home.html', task_id=task_id)





