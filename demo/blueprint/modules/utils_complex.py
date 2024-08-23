import os


def get_samples_list(select_raw_folder,fq_type,raw_myOptions):
    samples_list = []
    for filename in os.listdir(select_raw_folder):
        base_name = os.path.basename(filename)
        if filename.endswith(".fq.gz") and fq_type == "single" and raw_myOptions == "raw.fq.gz":
            samples_list.append(base_name[:-6])
        elif filename.endswith("_1.fq.gz") and fq_type == "pair" and raw_myOptions == "raw.fq.gz" and (
                filename[:-8] + "_2.fq.gz" in os.listdir(select_raw_folder)):
            samples_list.append(base_name[:-8])
        elif filename.endswith(".bam") and raw_myOptions == "raw.bam":
            samples_list.append(base_name[:-4])
    return samples_list


def find_rec_Option(type,options):
    if type=="fq1_fq2":
        paired_strings=[]
        fq1_strings = [s for s in options if s.endswith("_1.fq.gz")]
        fq2_strings = [s for s in options if s.endswith("_2.fq.gz")]
        for s1 in fq1_strings:
            s2 = s1[:-len("_1.fq.gz")] + "_2.fq.gz"
            if s2 in fq2_strings:
                paired_strings.append((s1, s2))
        if len(paired_strings) != 0 :
            return max(paired_strings, key=lambda pair: len(pair[0]) + len(pair[1]))
        else:
            return "", ""
    elif type=="fq1":
        fq1_strings = [s for s in options if s.endswith("_1.fq.gz")]
        if len(fq1_strings) != 0 :
            return max(fq1_strings, key=len)
        else:
            return ""
    elif type=="fq2":
        fq2_strings = [s for s in options if s.endswith("_2.fq.gz")]
        if len(fq2_strings) != 0 :
            return max(fq2_strings, key=len)
        else:
            return ""
    else:
        if len(options) != 0 :
            return max(options, key=len)
        else:
            return ""

def get_init_div(stepId, step_type, step_params, options_list, mode):
    options_fq = ''.join([f'<option value="{option}">{option}</option>' for option in options_list["fq"]])
    options_bam = ''.join([f'<option value="{option}">{option}</option>' for option in options_list["bam"]])
    options_other = ''.join([f'<option value="{option}">{option}</option>' for option in options_list["other"]])
    options_all = ''.join([f'<option value="{option}">{option}</option>' for option in options_list["all"]])
    
    
    start_sting = f"""
            <div id="input-step-{stepId}" class="form-group">
                <label><strong>step{stepId}: {step_type}:</strong></label>
                <input id="input-{stepId}" type="text" class="form-control" placeholder="输入额外参数（默认为空）" onchange="updateInput({stepId})" value="{step_params}">
            """
    end_string = f"""
                </div>
            """

    if mode == "single":
        if step_type == "hairpin_cut":
            result_sting = f"""
                <select id="select-{stepId}-fq1" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq1">{options_fq}</select>
                <select id="select-{stepId}-fq2" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq2">{options_fq}</select>
                <select id="select-{stepId}-fq" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq">{options_fq}</select>
            """
            rec_Option = {"fq1": find_rec_Option("fq1_fq2", options_list["fq"])[0], "fq2": find_rec_Option("fq1_fq2", options_list["fq"])[1], "fq": find_rec_Option("fq", options_list["fq"])}
        else:
            result_sting = f"""
                <select id="select-{stepId}-fq" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq">{options_fq}</select>
            """
            rec_Option = {"fq": find_rec_Option("fq", options_list["fq"])}

    else:
        if step_type == "hairpin_cut":
            result_sting = f"""
                <select id="select-{stepId}-fq1" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq1">{options_fq}</select>
                <select id="select-{stepId}-fq2" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq2">{options_fq}</select>
                <select id="select-{stepId}-fq" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq">{options_fq}</select>
            """
            rec_Option = {"fq1": find_rec_Option("fq1_fq2", options_list["fq"])[0], "fq2": find_rec_Option("fq1_fq2", options_list["fq"])[1], "fq": find_rec_Option("fq", options_list["fq"])}
        else:
            result_sting = f"""
                <select id="select-{stepId}-fq1" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq1">{options_fq}</select>
                <select id="select-{stepId}-fq2" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="fq2">{options_fq}</select>
            """
            rec_Option = {"fq1": find_rec_Option("fq1_fq2", options_list["fq"])[0], "fq2": find_rec_Option("fq1_fq2", options_list["fq"])[1]}
    return start_sting + result_sting + end_string , rec_Option

def get_result_option(step_id, step_type, select_option,my_options):
    sample_name="raw"
    if step_type.split()[0]=="Trim":
        temp="step"+step_id  + "_"+ step_type.split()[1] +"."+ sample_name + ".fq.gz"
    elif step_type.split()[0]=="map":
        temp="step"+step_id  + "_"+ step_type.split()[1] +"."+ sample_name + ".bam"
    elif step_type.split()[0]=="BAM":
        temp="step"+step_id  + "_"+ step_type.split()[1] +"."+ sample_name + ".bam"
    else:
        temp="step"+step_id  + "_"+ step_type.split()[1] +"."+ sample_name
    my_options.append(temp)
    return my_options


def get_type_index(step_type):
    if step_type.split()[1]=="bowtie2":
        return ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    elif step_type.split()[1]=="star":
        return ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    else:
        return ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
   



def get_type_option(my_options):
    fq=[]
    bam=[]
    other=[]
    for i in my_options:
        if i[-3:]==".fq" or i[-6:]==".fq.gz":
            fq.append(i)
        elif i[-3:]=="bam":
            bam.append(i)
        else:
            other.append(i)

    return {"fq":list(set(fq)),"bam":list(set(bam)),"other":list(set(other)),"all":list(set(my_options))}





tools_config = {
    "fastp": {
        "pair": (
            "      fq1 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz\"\n"
            "      fq2 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz\"\n" ),
        "single": "      fq1 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz\"\n",
    },
    "Trim": {
        "pair": (
            "      fq1 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz\"\n"
            "      fq2 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz\"\n"),
        "single": "      fq1 : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz\"\n",
    },
    "map": {
        "pair": "      bam : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam\"\n",
        "single": "      bam : \"{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam\"\n",
    },
}
# 这里的主要功能是 维护，每一步的输出，step_id反而是不需要的，前提是按照顺序
# 这里提供一个input的思路，先完成output的列表，再提取input
# 基于以上一点，这里要求input_list，必须开头是stepx_的形式，明确是哪一步的输出（原始输入除外），目前的批量版本，只需要raw开头验证原始即可


def generate_steps_config(soft_list,type_list,input_list, select_raw_folder, mode,raw_myOptions ):
    output_config_list=[]
    input_config_list = []
    for i in range(len(soft_list)):
        if soft_list[i] in tools_config.keys():
            output_config_list.append(tools_config[soft_list[i]][mode].format(input_step="step"+str(i),input_soft= soft_list[i] ))
        else:
            output_config_list.append(tools_config[type_list[i]][mode].format(input_step="step"+str(i), input_soft=soft_list[i]))
    for file in input_list:
        if file.startswith("raw"):
            if raw_myOptions[-3:] == "bam":
                input_config_list.append("      bam : \"" + select_raw_folder + "/{sample}.bam\"\n")
            elif mode=="single":
                input_config_list.append("      fq1 : \"" + select_raw_folder + "/{sample}.fq.gz\"\n")
            else :
                input_config_list.append(("      fq1 : \"" + select_raw_folder + "/{sample}_1.fq.gz\"\n"
                                          "      fq2 : \"" + select_raw_folder + "/{sample}_2.fq.gz\"\n"))
        else:
            input_config_list.append(output_config_list[int(file.split("_")[0][4:])])
    result= "steps:\n"
    for i in range(len(soft_list)):
        result += "  step" + str(i) + ":\n"
        result += "    input:\n"
        result += input_config_list[i]
        result += "    output:\n"
        result += output_config_list[i]

    return result, output_config_list




my_code_script = {
    "Trim fastp": {
        "pair": """
        mkdir -p $output_dir/${{sample}}/{step}_fastp
        $fastp {params} -i {filename}_1.fq.gz -I {filename}_2.fq.gz -o $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}_1.fq.gz \\
        -O $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}_2.fq.gz -j $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.json \\
        -h $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.html -w 4 >> $output_dir/logs/${{sample}}.log 2>&1
        """,
        "single": """
        mkdir -p $output_dir/${{sample}}/{step}_fastp
        $fastp {params} -i {filename}.fq.gz -o $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.fq.gz \\
        -j $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.json -h $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.html \\
        -w 4 >> $output_dir/logs/${{sample}}.log 2>&1
        """,
    },
    "Trim trim_galore": {
        "pair": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        trim_galore {params} --paired --cores 4 -o $output_dir/${{sample}}/{step}_trim_galore {filename}_1.fq.gz {filename}_2.fq.gz >> $output_dir/logs/${{sample}}.log 2>&1 
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {filename})_1_val_1.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}_1.fq.gz
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {filename})_2_val_2.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}_2.fq.gz  
        """,
        "single": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        trim_galore {params} --cores 4 -o $output_dir/${{sample}}/{step}_trim_galore {filename}.fq.gz >> $output_dir/logs/${{sample}}.log 2>&1 
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {filename})_trimmed.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}.fq.gz
        """
    }
}



# 这里要处理每一步的输入
# 情况1，原始输入    情况2，其他步骤的输入，构成为  $output_dir + sample + step +soft +  ，  重点是step+soft
# 这里对input_list要求变高
# 这里除了要steox之外，还要定位_raw
# 我决定修改input_list的规则，用stepx/soft.raw.fq.gz的形式


def generate_present_script(types_list,mode,params_list , input_list, select_raw_folder ):
    result=""
    for i in range(len(types_list)):
        if input_list[i].startswith("raw"):
            result += my_code_script[types_list[i]][mode].format(params=params_list[i], filename= select_raw_folder+"/${sample}", step="step"+str(i))
        else:
            prefix = str(input_list[i]).split("_")[0]
            suffix = input_list[i][len(prefix)+1:].split(".")[0]
            temp = f"$output_dir/${{sample}}/{prefix}_{suffix}/{suffix}_${{sample}}"
   
            result += my_code_script[types_list[i]][mode].format(params=params_list[i], filename=temp, step="step" + str(i))
    return result



available_steps = {
    'Trim 操作': ['Trim trim_galore', 'Trim fastp', 'Trim cutadapt'],
    '比对操作': ['map bowtie2', 'map bwa', 'map bismark'],
    'BAM 文件操作': ['BAM picard'],
    '其他': ['other qualimap', 'other featureCounts']
}

soft_format = {
    'trim_galore': "trim_galore  {params} --cores 4 -o ./ {input}.fq.gz",
    'fastp': "fastp {params} -i {input}.fq.gz -o {input}_trimmed.fq.gz -j {input}.json -h {input}.json -w 4",
}

