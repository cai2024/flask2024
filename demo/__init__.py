import os
import click
import json
from flask import Flask
from flask_bootstrap import Bootstrap4
from demo.blueprint.biology import biology
from demo.blueprint.analysis import analysis_bp
from demo.blueprint.complex_junction import complex_junction_bp

class BaseConfig:
    SECRET_KEY = os.getenv('SECRET_KEY', 'dev key')


def create_app(config_name=None):
    app = Flask('demo')
    app.config.from_object(BaseConfig)

    bootstrap=Bootstrap4()
    bootstrap.init_app(app)
    app.register_blueprint(biology)
    app.register_blueprint(analysis_bp)
    app.register_blueprint(complex_junction_bp)
    return app
