from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm







class WGBSForm_custom(MappingCheckForm):
    trim_soft = SelectField('选择trim软件',choices=["fastp","trim_galore"],validate_choice=False, validators=[DataRequired()],default="fastp")
    trim_params = StringField('trim额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    align_soft = SelectField('选择比对软件',choices=["bismark","bwa"],validate_choice=False, validators=[DataRequired()],default="bismark")
    align_params = StringField('比对额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="--non_directional")




    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='')



    dup_soft = SelectField('选择去重软件',choices=["picard","bwa"],validate_choice=False, validators=[DataRequired()],default="picard")
    dup_params =  StringField('去重额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")    
    call_modify_soft = SelectField('call_modify软件',choices=["bismark","bwa"],validate_choice=False, validators=[DataRequired()],default="bismark")
    call_modify_params = StringField('call_modeify额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    


  
    spikein_params = StringField('spikein额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')



class WGBSForm_default(MappingCheckForm):
    trim_soft = HiddenField(default="fastp")
    trim_params = HiddenField(default="-D")
    align_soft = HiddenField(default="bismark")
    align_params = HiddenField(default="--non_directional")


    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"],
                          validate_choice=False, validators=[DataRequired()], description='')


    dup_soft = HiddenField(default="picard")
    dup_params = HiddenField(default="")
    call_modify_soft = HiddenField(default="bismark")
    call_modify_params = HiddenField(default="")

    spikein_params = HiddenField(default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')




  
        








