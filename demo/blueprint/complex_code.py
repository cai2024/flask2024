from flask import request, Blueprint, render_template, current_app, url_for, redirect, Response, stream_with_context,g, session, Flask, jsonify,send_from_directory, abort
import os
import uuid
from datetime import datetime
from .modules.paths import paths_dict
complex_code = Blueprint('complex_code', __name__)




@complex_code.route('/complex_code/', methods=['GET', 'POST'])
def code_home():
    return render_template('complex_code/home.html', available_steps=available_steps)



available_steps = {
    'Trim 操作': ['Trim trim_galore', 'Trim fastp', 'Trim cutadapt'],
    '比对操作': ['map bowtie2', 'map bwa', 'map bismark'],
    'BAM 文件操作': ['BAM picard'],
    '其他': ['other qualimap', 'other featureCounts']
}


soft_format={
    'trim_galore': "trim_galore  {params} --cores 4 -o ./ {input}.fq.gz",
    'fastp': "fastp {params} -i {input}.fq.gz -o {input}_trimmed.fq.gz -j {input}.json -h {input}.json -w 4",
    'cutadapt': "cutadapt -o {params} cutadapt_{input}.fq.gz {input}.fq.gz",
    'bowtie2': "bowtie2 {params} -x /index -U {input}.fq.gz -S {input}.sam  -p 10 && samtools view -b {input}.sam | samtools sort -@ 10 -o {input}.bam",
    'bwa': "bwa mem {params} -t 10 /index {input}.fq.gz | samtools sort -@ 10 -o {input}.bam",
    'bismark': "bismark {params} --path_to_bowtie2 /bowtie2 /index {input}.fq.gz -o ./ -p 10 && samtools sort -@ 10  bismark.bam > {input}.bam",
    'picard': "java -jar picard MarkDuplicates --REMOVE_DUPLICATES {params} -I {input}.bam -O de_{input}.bam -M {input}.report",
    'qualimap': "qualimap",
    'featureCounts': "featureCounts",
        }





soft_script_single={
    "Trim trim_galore": 
    """
$trim_galore  {params} --cores 4 -o $output_dir/{step}_trim_galore  ${fastq1}.fq.gz  
mv $output_dir/{step}_trim_galore/${fastq1}_trimmed.fq.gz $output_dir/{step}/trim_galore_${fastq1}.fq.gz
    """,

    "Trim fastp":
    """
$trim_galore  {params} --cores 4 -o $output_dir/{step}_trim_galore  ${fastq1}.fq.gz  
mv $output_dir/{step}_trim_galore/${fastq1}_trimmed.fq.gz $output_dir/{step}/trim_galore_${fastq1}.fq.gz
    """,

    "Trim cutadapt":
    """
$trim_galore  {params} --cores 4 -o $output_dir/{step}_trim_galore  ${fastq1}.fq.gz  
mv $output_dir/{step}_trim_galore/${fastq1}_trimmed.fq.gz $output_dir/{step}/trim_galore_${fastq1}.fq.gz
    """,

    'map bowtie2': 
    """
bowtie2 {params} -x /index -U {fastq1}.fq.gz -S {fastq1}.sam  -p 10 
samtools view -b {fastq1}.sam | samtools sort -@ 10 -o {fastq1}.bam
    """,

    'map bwa':
    """
bowtie2 {params} -x /index -U {fastq1}.fq.gz -S {fastq1}.sam  -p 10 
samtools view -b {fastq1}.sam | samtools sort -@ 10 -o {fastq1}.bam
    """,

    'map bismark':
    """
bowtie2 {params} -x /index -U {fastq1}.fq.gz -S {fastq1}.sam  -p 10 
samtools view -b {fastq1}.sam | samtools sort -@ 10 -o {fastq1}.bam
    """,

    'BAM picard': 
    """
"picard"
    """,
}










@complex_code.route('/get_preset', methods=['POST'])
def preset_soft():
    data = request.json
    presetType=data.get("presetType")
    if presetType=="WGS":
        options=['Trim trim_galore','map bowtie2', 'BAM picard']
        select=['raw.fq','step0_raw.fq', 'step1_step0_raw.fq']
        params=["-A NNNN","","-D"]

    return jsonify({"options":options,"select":select,"params":params})



def get_type_option(step_type,my_options):
    fq=[]
    bam=[]
    other=[]
    for i in my_options:
        if i[-3:]==".fq":
            fq.append(i)
        elif i[-3:]=="bam":
            bam.append(i)
        else:
            other.append(i)
    if step_type.split()[0]=="Trim":
        return list(set(fq))
    elif step_type.split()[0]=="map":
        return list(set(fq))
    elif step_type.split()[0]=="BAM":
        return list(set(bam))
    else:
        return list(set(other))

def get_result_option(step_id, step_type, select_option,my_options):
    if step_type.split()[0]=="Trim":
        temp="step"+step_id  + "_"+ select_option
    elif step_type.split()[0]=="map":
        temp="step"+step_id  + "_"+ select_option[:-2]+"bam"
    elif step_type.split()[0]=="BAM":
        temp="step"+step_id  + "_"+ select_option[:-3]+"bam"
    else:
        temp="step"+step_id  + "_"+ select_option
    my_options.append(temp)
    return my_options
        




@complex_code.route('/submit_steps',  methods=['POST'])
def submit_steps():
    task_id = f"{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())}"
    data = request.json
    steps = data.get('steps')
    fq_type=data.get('fq_type')
    raw_file=data.get('raw_file')

    head_text=f"""
sample="{raw_file.split('.')[0]}"
mode="{fq_type}"
output="/data/output/{task_id}"

"""
    print(fq_type)
    print(raw_file)
    result_list=[]

    all_soft=[]
    for dic in steps:
        step_id = str(dic["id"])
        step_type = dic["type"]
        step_params = dic["oninput"]
        step_file = dic["select_option"]
        all_soft.append(f'{step_type.split()[1]}="{paths_dict[step_type.split()[1]]}"') 

        temp=soft_script_single[step_type].format(params=step_params,fastq1=step_file,step=step_id)
        result_list.append(temp)
    
    final_list = all_soft + result_list

    return jsonify({'my_code': head_text+'\n'.join(final_list),"task_id":task_id})
   
@complex_code.route('/submit_form_code/<task_id>')
def submit_form_code(task_id):
    return render_template('complex_code/code.html')





@complex_code.route('/get_options', methods=['POST'])
def get_options():
    data = request.json
    step_type = data.get('step_type')
    my_options = data.get('my_options')
    options = get_type_option(step_type,my_options)
    return jsonify({'options': options})



@complex_code.route('/update_my_options', methods=['POST'])
def update_my_options():
    data = request.json
    step_id = str(data.get('step_id'))
    step_type = data.get('step_type')
    select_option = data.get('select_option')
    my_options = data.get('my_options')
    
    
    my_options=get_result_option(step_id,step_type,select_option,my_options)
    my_options=list(set(my_options))

    return jsonify({'my_options': my_options})


@complex_code.route('/get_display_value', methods=['POST'])
def update_display_params():
    data = request.get_json()
    step_type=data['type'].split()[1]
    step_params=data['params']
    select_option = data['select_option']
    
    if select_option[-3:]==".fq":
        input_file=select_option[:-3]
    else:
        input_file=select_option[:-4]
    displayString = soft_format[step_type].format(params=step_params,input=input_file)

    # 假设这个处理结果就是你要返回的展示字符串
    return jsonify({'displayString': displayString})

@complex_code.route('/index_select', methods=['POST'])
def update_index():
    data = request.json
    step_type = data.get('step_type').split()[1]
    index_options=["mouse mm9", "mouse mm10", "mouse mm39","human hg18","human hg19","human hg38","human hs1"]
    return jsonify({'index_options': index_options})
        



