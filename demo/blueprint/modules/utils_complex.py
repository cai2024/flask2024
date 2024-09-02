import os

available_steps = {
    'Trim 操作': ['Trim trim_galore', 'Trim fastp', 'Trim cutadapt'],
    '比对操作': ['map bowtie2', 'map bwa', 'map bismark'],
    'BAM 文件操作': ['BAM picard'],
    '其他': ['other qualimap', 'other featureCounts']
}
my_code_script = {
    "Trim fastp": {
        "pair": """
        mkdir -p $output_dir/${{sample}}/{step}_fastp
        $fastp {params} -i {fq1} -I {fq2} -o $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}_1.fq.gz \\
        -O $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}_2.fq.gz -j $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.json \\
        -h $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.html -w 4 >> $output_dir/logs/${{sample}}.log 2>&1
        """,
        "single": """
        mkdir -p $output_dir/${{sample}}/{step}_fastp
        $fastp {params} -i {fq} -o $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.fq.gz \\
        -j $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.json -h $output_dir/${{sample}}/{step}_fastp/fastp_${{sample}}.html \\
        -w 4 >> $output_dir/logs/${{sample}}.log 2>&1
        """,
    },
    "Trim trim_galore": {
        "pair": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        trim_galore {params} --paired --cores 4 -o $output_dir/${{sample}}/{step}_trim_galore {fq1} {fq2} >> $output_dir/logs/${{sample}}.log 2>&1 
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {fq1} _1.fq.gz)_1_val_1.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}_1.fq.gz
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {fq2} _2.fq.gz)_2_val_2.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}_2.fq.gz  
        """,
        "single": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        trim_galore {params} --cores 4 -o $output_dir/${{sample}}/{step}_trim_galore {fq} >> $output_dir/logs/${{sample}}.log 2>&1 
        mv $output_dir/${{sample}}/{step}_trim_galore/$(basename {fq} .fq.gz)_trimmed.fq.gz $output_dir/${{sample}}/{step}_trim_galore/trim_galore_${{sample}}.fq.gz
        """
    }
}
my_code_simple = {
    "Trim fastp": {
        "pair": "fastp {params} -i sample_1.fq.gz -I sample_2.fq.gz -o fastp_sample_1.fq.gz -O fastp_sample_2.fq.gz -j fastp_sample.json -h fastp_sample.html -w 4 ",
        "single": "fastp {params} -i sample.fq.gz -o fastp_sample.fq.gz -j fastp_sample.json -h fastp_sample.html -w 4 ",
    },
    "Trim trim_galore": {
        "pair": "trim_galore {params} --paired --cores 4 -o step{stepId}_trim_galore sample_1.fq.gz sample_2.fq.gz ",
        "single": "trim_galore {params} --cores 4 -o step{stepId}_trim_galore sample.fq.gz "
    }
}


tool_input = {
    "fastp": {
        "pair": {"fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" ,
                 "fq2": "step{step_id}_{step_type}.{sample}_2.fq.gz" },
        "single": {"fq": "step{step_id}_{step_type}.{sample}.fq.gz"}
    },
    "Trim": {
        "pair": {"fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" ,
                 "fq2": "step{step_id}_{step_type}.{sample}_2.fq.gz" },
        "single": {"fq": "step{step_id}_{step_type}.{sample}.fq.gz"}
    },
    "map": {
        "pair": {"fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" ,
                 "fq2": "step{step_id}_{step_type}.{sample}_2.fq.gz" },
        "single": {"fq": "step{step_id}_{step_type}.{sample}.fq.gz"}
    },
}
tools_config = {
    "fastp": {
        "pair": {"fq1": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz",
                 "fq2": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz"},
        "single": {"fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"}
    },
    "Trim": {
        "pair": {"fq1": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz",
                 "fq2": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz"},
        "single": {"fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"}
    },
    "map": {
        "pair": {"bam": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam"},
        "single": {"bam": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam"}
    },
}

def get_init_div(stepId, step_type, step_params, options_list, mode):
    option_dir = { key: ''.join(f'<option value="{option}">{option}</option>' for option in options)
        for key, options in options_list.items() }
    result_string = f"""
            <div id="input-step-{stepId}" class="form-group">
                <label><strong>step{stepId}: {step_type}:</strong></label>
                <input id="input-{stepId}" type="text" class="form-control" placeholder="输入额外参数（默认为空）" onchange="updateInput({stepId})" oninput="updateOutputArea()" value="{step_params}">
            """
    key_list = tools_config.get(step_type.split()[1], tools_config[step_type.split()[0]])[mode].keys()
    result_string += "\n".join(
        f'<select id="select-{stepId}-{title}" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="{title}">{option_dir.get(title, [])}</select>'
        for title in key_list)
    result_string += f"\n </div> \n"
    key_dir = {key :options_list[key][max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                for temp in options_list[key]]),key=lambda x: x[1])[0]] for key in key_list }
    return  result_string , key_dir

def get_result_option(step_id, step_type, mode, temp_options):
    key_list = tool_input.get(step_type.split()[1], tool_input[step_type.split()[0]])[mode].keys()
    for key in key_list:
        temp_options[key].append(tool_input.get(step_type.split()[1], tool_input[step_type.split()[0]])[mode][key]
                                 .format(step_id=step_id,step_type=step_type.split()[1],sample="raw"))
    temp_options["all"].extend(tool_input.get(step_type.split()[1], tool_input[step_type.split()[0]])[mode].values())
    return temp_options

def get_type_index(stepId, step_type):
    if step_type.split()[1] == "bowtie2":
        select_list = ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    elif step_type.split()[1] == "star":
        select_list = ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    else:
        select_list = ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    optionsHtml = ''.join([f'<option value="{option}">{option}</option>' for option in select_list])
    return f"""
            <div id="index-step-${stepId}" class="form-group">
                <label><strong>step${stepId}: ${step_type} 参考基因组:</strong></label>
                <select class="form-control mt-2">${optionsHtml}</select>
            </div>       
            """


def update_select_dir(key_dir, change_key, step_type, mode):
    if step_type == "hairpin_cut" and mode == "single":
        pass
    elif change_key == "fq1" and mode == "pair" and ("fq2" in key_dir):
        paired_file = key_dir["fq1"][0].replace('_1.fq.gz', '_2.fq.gz')
        key_dir["fq2"][0] = paired_file if paired_file in key_dir["fq2"][1] else key_dir["fq2"][0]
    elif change_key == "fq2" and mode == "pair" and ("fq1" in key_dir):
        paired_file = key_dir["fq2"][0].replace('_2.fq.gz', '_1.fq.gz')
        key_dir["fq1"][0] = paired_file if paired_file in key_dir["fq1"][1] else key_dir["fq1"][0]
    return key_dir




def generate_steps_config(steps, select_raw_folder, mode, task_id, samples_list):
    input_list = [record["select_option"] for record in steps]  # 这是字典
    soft_list = [record["type"].split()[1] for record in steps]
    type_list = [record["type"].split()[0] for record in steps]
    params_list = [record["params"] for record in steps]

    output_config_list = []
    for i, soft in enumerate(soft_list):
        config = tools_config.get(soft, tools_config[type_list[i]])[mode]
        output_config_list.append({k: v.format(input_step=f"step{i}", input_soft=soft) for k, v in config.items()})

    input_config_list = []
    for file_dir in input_list:
        temp_dir = {}
        for k, v in file_dir.items():
            if v.startswith("raw"):
                suffix = {"fq1": "_1.fq.gz", "fq2": "_2.fq.gz", "fq": ".fq.gz"}.get(k, f".{k}")
                temp_dir[k] = f"{select_raw_folder}/{{sample}}{suffix}"
            else:
                temp_dir[k] = output_config_list[int(v.split("_")[0][4:])][k]
        input_config_list.append(temp_dir)

    config_steps = "steps:\n"
    for i in range(len(soft_list)):
        config_steps += "  step" + str(i) + ":\n"
        config_steps += "    input:\n"
        config_steps += "\n".join(f"      {key} : \"{value}\"" for key, value in input_config_list[i].items()) + "\n"
        config_steps += "    output:\n"
        config_steps += "\n".join(f"      {key} : \"{value}\"" for key, value in output_config_list[i].items()) + "\n"

    config_string = (f"mode: {mode}\n"
                     f"output_dir: /data/cailab/flask_out/{task_id}\n")
    config_string += "samples:\n" + "\n".join(f"  - {sample}" for sample in samples_list) + "\n"
    config_string += "\n".join(f"step{i}_params : \"{params}\"" for i, params in enumerate(params_list)) + "\n"
    config_string += config_steps

    snakemake_string = (
        "configfile: \"config.yaml\"\n"
        "include: \"/data/cailab/flask2024/smk/trim.smk\"\n"
        "rule all:\n"
        "    input:\n")
    snakemake_string += "\n".join(
        f"        expand(\"{temp}\", sample=config[\"samples\"], output_dir=config[\"output_dir\"])," for temp in
        [list(record.values())[0] for record in output_config_list]) + "\n"

    return config_string, snakemake_string


def generate_present_script(steps, mode, select_raw_folder, samples_list, task_id, paths_dict):
    input_list = [record["select_option"] for record in steps]  # 这是字典
    soft_list = [record["type"].split()[1] for record in steps]
    types_list = [record["type"] for record in steps]
    params_list = [record["params"] for record in steps]

    formatted_string = ""
    for i, type in enumerate(types_list):
        values_dict = {}
        for key, value in input_list[i].items():
            if value.startswith("raw"):
                values_dict[key] = select_raw_folder + "/${sample}" + value.split("raw", 1)[1]
            else:
                prefix = str(value).split("_")[0]
                suffix = value[len(prefix) + 1:].split(".")[0]
                values_dict[key] = f"$output_dir/${{sample}}/{prefix}_{suffix}/{suffix}_${{sample}}" + \
                                   value.split(".raw", 1)[1]

        values_dict["params"] = params_list[i]
        values_dict["step"] = "step" + str(i)
        formatted_string += my_code_script[type][mode].format(**values_dict)

    my_code = (f"mode=\"{mode}\"\n"
               f"output_dir=\"/data/output/{task_id}\"\n")
    my_code += '\n'.join([f'{soft}="{paths_dict[soft]}"' for soft in soft_list]) + "\n"
    my_code += 'samples=("{}")'.format('" "'.join(samples_list)) + "\n"
    my_code += 'mkdir -p $output_dir/logs\n'
    my_code += 'for sample in "${samples[@]}"; do\n'
    my_code += formatted_string + "\n"
    my_code += "done\n"

    return my_code






