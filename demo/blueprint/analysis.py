from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os

analysis_bp = Blueprint('analysis', __name__)



@analysis_bp.route('/docs/<path:path>')
def serve_doc(path):
    base_directory = 'templates/analysis'
    if path.endswith('/'):
        # 请求目录的 index.html
        return send_from_directory(base_directory, f"{path}/index.html")
    else:
        # 直接返回文件
        return send_from_directory(base_directory, path)

@analysis_bp.route('/docs/')
def index():
    # 默认服务 index.html
    return send_from_directory('templates/analysis', 'index.html')



@analysis_bp.route('/analysis')
def analysis():
    return redirect(url_for('analysis.serve_doc', path=''))
