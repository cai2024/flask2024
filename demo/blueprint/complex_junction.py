from flask import Blueprint, render_template

complex_junction_bp = Blueprint('complex_junction', __name__, template_folder='templates')

@complex_junction_bp.route('/complex_junction')
def complex_junction():
    return render_template('biology/mapping.html')

