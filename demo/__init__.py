import os
import click
import json
from flask import Flask
from flask_bootstrap import Bootstrap4
from demo.blueprint.biology import biology
from demo.blueprint.analysis import analysis_bp
from demo.blueprint.hidden import hidden_bp
from demo.blueprint.batch_code import batch_code
from demo.blueprint.individual_code import individual_code
from demo.blueprint.common import common_blueprint
from demo.blueprint.modules.utils_common import load_config


def create_app(config_name=None):
    app = Flask('demo')
    app.config['SECRET_KEY'] = 'dev key'
    app.config['DEBUG'] = True
    app.config['DATABASE_URI'] = 'sqlite:///mydatabase.db'

    available_steps, soft_config, index_options_dict, paths_dict, code_script = load_config()
    app.config['available_steps'] = available_steps
    app.config['soft_config'] = soft_config
    app.config['index_options_dict'] = index_options_dict
    app.config['paths_dict'] = paths_dict
    app.config['code_script'] = code_script

    bootstrap=Bootstrap4()
    bootstrap.init_app(app)
    app.register_blueprint(biology)
    app.register_blueprint(hidden_bp)
    app.register_blueprint(analysis_bp)
    app.register_blueprint(batch_code)
    app.register_blueprint(individual_code)
    app.register_blueprint(common_blueprint)
    return app
