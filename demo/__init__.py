import os
import click
import json
from flask import Flask
from flask_bootstrap import Bootstrap4
from demo.blueprint.biology import biology
from demo.blueprint.analysis import analysis_bp
from demo.blueprint.hidden import hidden_bp
from demo.blueprint.batch_code import complex_code
from demo.blueprint.individual_code import individual_code
from demo.blueprint.common import common_blueprint

class BaseConfig:
    SECRET_KEY = os.getenv('SECRET_KEY', 'dev key')


def create_app(config_name=None):
    app = Flask('demo')
    app.config.from_object(BaseConfig)

    bootstrap=Bootstrap4()
    bootstrap.init_app(app)
    app.register_blueprint(biology)
    app.register_blueprint(hidden_bp)
    app.register_blueprint(analysis_bp)
    app.register_blueprint(complex_code)
    app.register_blueprint(individual_code)
    app.register_blueprint(common_blueprint)
    return app
