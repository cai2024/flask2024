from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm





class WGBSForm_default(MappingCheckForm):
    ref_ana = StringField('选择参考基因组物种和版本')
    submit_ana = SubmitField('开始比对')


class WGBSForm_custom(MappingCheckForm):
    ref_ana = StringField('选择参考基因组物种和版本')
    submit_ana = SubmitField('开始比对')
