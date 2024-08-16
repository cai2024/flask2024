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



def get_type_option(step_type,my_options):
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
    if step_type.split()[0]=="Trim":
        return list(set(fq))
    elif step_type.split()[0]=="map":
        return list(set(fq))
    elif step_type.split()[0]=="BAM":
        return list(set(bam))
    else:
        return list(set(other))



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

