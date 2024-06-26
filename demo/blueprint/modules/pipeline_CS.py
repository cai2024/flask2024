from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm







class CSForm_custom(MappingCheckForm):
    trim_soft = SelectField('选择trim软件',choices=["fastp","trim_galore"],validate_choice=False, validators=[DataRequired()],default="fastp")
    trim_params = StringField('trim额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    cut_soft = SelectField('选择hairpin cut软件',choices=["CGS_haipin","five or six"],validate_choice=False, validators=[DataRequired()],default="fastp")
    cut_params = StringField('hairpin cut额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")    
    align_soft = SelectField('选择比对软件',choices=["bowtie2","bwa"],validate_choice=False, validators=[DataRequired()],default="bowtie2")
    align_params = StringField('bowtie2额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")



    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='') 

    dup_soft = SelectField('选择去重软件',choices=["picard","samtools"],validate_choice=False, validators=[DataRequired()],default="picard")
    dup_params = StringField('去重参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    call_modify_soft = SelectField('call_modify软件',choices=["CGS_call_modify","six_call_modify"],validate_choice=False, validators=[DataRequired()],default="CGS_call_modify")
    call_modify_params = StringField('call_modeify额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")


    spikein_params = StringField('spikein额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')



class CSForm_default(MappingCheckForm):
    trim_soft = HiddenField(default="fastp")
    trim_params = HiddenField(default="-D")
    cut_soft = HiddenField(default="CGS_haipin")
    cut_params = HiddenField(default="")
    align_soft = HiddenField(default="bowtie2")
    align_params = HiddenField(default="")


    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='')


    dup_soft = HiddenField(default="picard")
    dup_params = HiddenField(default="")
    call_modify_soft = HiddenField(default="CGS_call_modify")
    call_modify_params = HiddenField(default="")

    spikein_params = HiddenField( default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')




  
        








