from flask import Blueprint, render_template

analysis_bp = Blueprint('analysis', __name__, template_folder='templates')

@analysis_bp.route('/analysis')
def analysis():
    return render_template('analysis/code.html')

