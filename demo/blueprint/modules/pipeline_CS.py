from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm







class CSForm_custom(MappingCheckForm):
    trim_soft = SelectField('选择trim软件',choices=["fastp","trim_galore"],validate_choice=False, validators=[DataRequired()],default="fastp")
    trim_params = StringField('trim额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    align_soft = SelectField('选择比对软件',choices=["bowtie2","bwa"],validate_choice=False, validators=[DataRequired()],default="bowtie2")
    align_params = StringField('bowtie2额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")



    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","huamn hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='') 

    hairpin_cut_params = StringField('picard额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="--min 50 --rule 2")
    picard_params = StringField('picard额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="--MINIMUM_DISTANCE -1")
    spikein_params = StringField('spikein额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')



class CSForm_default(MappingCheckForm):
    trim_soft = HiddenField(default="fastp")
    trim_params = HiddenField(default="-D")
    align_soft = HiddenField(default="bowtie2")
    align_params = HiddenField(default="")



    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","huamn hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='')
    hairpin_cut_params = HiddenField( default="--min 50 --rule 2")
    picard_params = HiddenField(default="")
    spikein_params = HiddenField( default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')




  
        








