from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm







class CSForm_custom(MappingCheckForm):
    fastp_params = StringField('fastp额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="-D")
    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","huamn hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='') 
    bowtie2_params = StringField('bowtie2额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")   
    picard_params = StringField('picard额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="--MINIMUM_DISTANCE 200")
    hairpin_cut_params = StringField('picard额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="--min 50 --rule 2") 
    spikein_params = StringField('spikein额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')



class CSForm_default(MappingCheckForm):
    fastp_params = HiddenField(default="-D")
    bowtie2_params = HiddenField(default="")
    picard_params = HiddenField(default="")
    hairpin_cut_params = HiddenField( default="--min 50 --rule 2")
    spikein_params = HiddenField( default="lambda,mc,G5hmc,puc19")
    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","huamn hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='')
    submit_ana = SubmitField('开始比对')




  
        








