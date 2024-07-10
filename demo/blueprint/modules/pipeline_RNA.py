from wtforms.validators import DataRequired, Length, Optional, URL
from wtforms import SubmitField, SelectField, BooleanField, StringField, IntegerField, RadioField, HiddenField
from flask_wtf import FlaskForm
from .utils import MappingCheckForm







class RNAForm_custom(MappingCheckForm):
    trim_soft = SelectField('选择trim软件',choices=["fastp","trim_galore"],validate_choice=False, validators=[DataRequired()],default="trim_galore")
    trim_params = StringField('trim额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
    align_soft = SelectField('选择比对软件',choices=["bowtie2","bwa","star"],validate_choice=False, validators=[DataRequired()],default="star")
    align_params = StringField('比对软件额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="")


    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=["mouse mm10", "mouse mm39","human hg19","human hg38"],
                          validate_choice=False, validators=[DataRequired()], description='') 
    
    ref_len = SelectField('选择索引长度(readlength-1)',choices=[ "74","99","149"],validate_choice=False, validators=[DataRequired()], description='')

    count_soft = SelectField('选择count软件',choices=["featureCounts"],validate_choice=False, validators=[DataRequired()],default="featureCounts")
    count_params = StringField('参数设置', validators=[Optional(), Length(min=0, max=50)], default="")
  
    spikein_params = StringField('spikein额外参数设置', validators=[Optional(), Length(min=0, max=50)], default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')



class RNAForm_default(MappingCheckForm):
    trim_soft = HiddenField(default="trim_galore")
    trim_params = HiddenField(default="")
    align_soft = HiddenField(default="star")
    align_params = HiddenField(default="")


    ref_ana = SelectField('选择参考基因组物种和版本',
                          choices=[ "mouse mm10", "mouse mm39","human hg19","human hg38"],
                          validate_choice=False, validators=[DataRequired()], description='')

    ref_len = SelectField('选择索引长度(readlength-1)',choices=[ "74","99","149"],validate_choice=False, validators=[DataRequired()], description='')

    count_soft = HiddenField(default="featureCounts")
    count_params = HiddenField(default="")

    spikein_params = HiddenField(default="lambda,mc,G5hmc,puc19")
    submit_ana = SubmitField('开始比对')




  
        








