import os

available_steps = {
    'Trim 操作': ['Trim trim_galore', 'Trim fastp', 'Trim cutadapt' , 'Trim hairpin_cut', 'Trim hairpin_cut_N'],
    '比对操作': ['map bowtie2', 'map bowtie2_spikein', 'map bwa','map bwa_spikein', 'map bismark','map bismark_spikein', 'map STAR' ],
    'BAM 文件操作': ['BAM remove_duplicates','BAM CGS_call_modify','BAM CGS_spikein_modify','BAM TBS_call_modify','BAM TBS_spikein_modify'],
    '其他': ['other qualimap_cover', 'other featureCounts']
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
    },
    "map bowtie2": {
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
    },
    "map bwa": {
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
    },
    "map bismark": {
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
    },
    "map STAR": {
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
    },
    "BAM remove_duplicates": {
        "pair": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        """,
        "single": """
        mkdir -p $output_dir/${{sample}}/{step}_trim_galore
        """
    },




}
my_code_simple = {
    "Trim fastp": {
        "pair": "fastp {params} -i sample_1.fq.gz -I sample_2.fq.gz -o fastp_sample_1.fq.gz -O fastp_sample_2.fq.gz -j fastp_sample.json -h fastp_sample.html -w 4 ",
        "single": "fastp {params} -i sample.fq.gz -o fastp_sample.fq.gz -j fastp_sample.json -h fastp_sample.html -w 4 ",
    },
    "Trim trim_galore": {
        "pair": "trim_galore {params} --paired --cores 4 -o step{stepId}_trim_galore sample_1.fq.gz sample_2.fq.gz ",
        "single": "trim_galore {params} --cores 4 -o step{stepId}_trim_galore sample.fq.gz "
    },
    "map bowtie2": {
        "pair": "trim_galore {params} --paired --cores 4 -o step{stepId}_trim_galore sample_1.fq.gz sample_2.fq.gz ",
        "single": "trim_galore {params} --cores 4 -o step{stepId}_trim_galore sample.fq.gz "
    },

}


tools_input = {
    "CGS_call_modify": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam" ,
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" },
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam",
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz"}
    },
    "CGS_spikein_modify": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam" ,
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" },
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam",
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz"}
    },
    "TBS_call_modify": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam" ,
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" },
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam",
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz"}
    },
    "TBS_spikein_modify": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam" ,
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz" },
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam",
                 "fq1": "step{step_id}_{step_type}.{sample}_1.fq.gz"}
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
    "BAM": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam"},
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam"}
    },
    "qualimap_cover": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam"},
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam"}
    },
    "other": {
        "pair": {"bam": "step{step_id}_{step_type}.{sample}.bam"},
        "single": {"bam": "step{step_id}_{step_type}.{sample}.bam"}
    },

}
tools_config = {
    "hairpin_cut": {
        "pair": {"fq1": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz",
                 "fq2": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz",
                  "fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"},
        "single": {"fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"}
    },
    "hairpin_cut_N": {
        "pair": {"fq1": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_1.fq.gz",
                 "fq2": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}_2.fq.gz",
                  "fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"},
        "single": {"fq": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.fq.gz"}
    },

    "CGS_call_modify": {
        "pair": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                 "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed",},
        "single": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                   "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed"}
    },
    "CGS_spikein_modify": {
        "pair": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                 "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed",},
        "single": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                   "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed"}
    },

    "TBS_call_modify": {
        "pair": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                 "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed",},
        "single": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                   "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed"}
    },
    "TBS_spikein_modify": {
        "pair": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                 "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed",},
        "single": {"bed_mc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.mc.all.bed",
                   "bed_hmc": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.hmc.all.bed"}
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
    "BAM": {
        "pair": {"bam": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam"},
        "single": {"bam": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.bam"}
    },
    "qualimap_cover": {
        "pair": {"other": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.txt"},
        "single": {"other": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.txt"}
    },
    "other": {
        "pair": {"other": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.txt"},
        "single": {"other": "{{output_dir}}/{{sample}}/{input_step}_{input_soft}/{input_soft}_{{sample}}.txt"}
    },


}

index_options_dict = {
    "bowtie2": {
        "ref": ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    },
    "STAR": {
        "ref": ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"],
        "gtf": ["a1", "a2"],
        "overhang": ["149", "99"]
    },
    "default": {
        "ref": ["mouse mm9", "mouse mm10", "mouse mm39", "human hg18", "human hg19", "human hg38", "human hs1"]
    }
}


def get_init_div(stepId, step_type, step_params, options_list, mode):
    # 每项软件的初始化包含的选项数量和类型，由tools_input.keys提供
    # 每个类型提供的选项，由options_list提供（字典）
    # 推荐选项获取step后最大值
    # 测试验证时间：2024.09.03
    option_dir = { key: ''.join(f'<option value="{option}">{option}</option>' for option in options) for key, options in options_list.items()}
    result_string = f"""
            <div id="input-step-{stepId}" class="form-group">
                <label><strong>step{stepId}: {step_type}:</strong></label>
                <input id="input-{stepId}" type="text" class="form-control" placeholder="输入额外参数（默认为空）" onchange="updateInput({stepId})" oninput="updateOutputArea()" value="{step_params}">
            """
    key_list = tools_input.get(step_type.split()[1], tools_input[step_type.split()[0]])[mode].keys()
    result_string += "\n".join(
        f'<select id="select-{stepId}-{title}" class="form-control mt-2" onchange="updateSelectOption({stepId})" title="{title}">{option_dir.get(title, [])}</select>'
        for title in key_list)
    result_string += f"\n </div> \n"
    key_dir = {key :options_list[key][max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                for temp in options_list[key]]),key=lambda x: x[1])[0]] for key in key_list }
    return  result_string , key_dir

def get_result_option(step_id, step_type, mode, temp_options):
    # 每次扩充当前id的软件的生成文件
    # 逻辑上是先获取目标软件生成的文件类型有哪些,通过suffix模式扩充temp_options字典
    # 测试验证时间：2024.09.02
    # 原始逻辑没有考虑同一个软件生成多个同类型文件，对文件名后缀的要求太死板
   
    key_list = list(tools_config.get(step_type.split()[1], tools_config[step_type.split()[0]])[mode].keys())
    for key in key_list:
        temp_string = tools_config.get(step_type.split()[1], tools_config[step_type.split()[0]])[mode][key].split("/")[-1].replace("_{{sample}}","_raw")
        temp_options[key.split("_")[0]].append(f"step{step_id}_" + temp_string.format(input_soft=step_type.split()[1]))
        temp_options["all"].append(f"step{step_id}_" + temp_string.format(input_soft=step_type.split()[1]))
    return temp_options

def get_type_index(stepId, step_type):
    # index定义为字典形式，id为index-stepId-key
    # 测试验证时间：2024.09.03
    select_dir = index_options_dict.get(step_type.split()[1], index_options_dict["default"])
    result= f'<div id="index-step-{stepId}" class="form-group">\n'
    for title in select_dir:
        option_div = ''.join([f'<option value="{option}">{option}</option>' for option in select_dir.get(title, [])])
        result += f'<select id="index-{stepId}-{title}" class="form-control mt-2" title="{title}">{option_div}</select>\n'
    return f"<label><strong>step{stepId}: {step_type} 参考基因组:</strong></label>\n{result}\n</div>\n"

def update_select_dir(key_dir, change_key, step_type, mode):
    # key_dir[key]=[step.select_option[key],optionsList],step.select_option[key]为当前选项，optionsList为可选选项
    # 逻辑是在js创建上述格式字典，修改change_key配对的key的step.select_option[key]
    # 返回后js将选项和step做出对应修改
    # 测试验证时间：2024.09.03
    if step_type == "hairpin_cut" and mode == "single":
        pass
    elif change_key == "fq1" and mode == "pair" and ("fq2" in key_dir):
        paired_file = key_dir["fq1"][0].replace('_1.fq.gz', '_2.fq.gz')
        key_dir["fq2"][0] = paired_file if paired_file in key_dir["fq2"][1] else key_dir["fq2"][0]
    elif change_key == "fq2" and mode == "pair" and ("fq1" in key_dir):
        paired_file = key_dir["fq2"][0].replace('_2.fq.gz', '_1.fq.gz')
        key_dir["fq1"][0] = paired_file if paired_file in key_dir["fq1"][1] else key_dir["fq1"][0]
    return key_dir

def generate_steps_config(steps, select_raw_folder, mode, task_id, samples_list, paths_dict):
    # 对于config，1.index读入config，包含ref_species，ref_release和other
    # 2.mode，output_dir，soft路径（待添加）
    # 3.samples 4.params 5.steps
    # 对于snakemake，1.configfile 2.include全部的smk 3.自动获取输出
    index_list = [record["select_index"] for record in steps if  record["select_index"]]
    input_list = [record["select_option"] for record in steps]  
    soft_list = [record["type"].split()[1] for record in steps]
    type_list = [record["type"].split()[0] for record in steps]
    params_list = [record["params"] for record in steps]
    mode_list = [record["mode"] for record in steps]


    config_dir={}
    # config_dir["mode"]=mode
    config_dir["output_dir"]=f"/data/cailab/flask_out/{task_id}"
    config_dir.update(paths_dict)
    for temp in index_list:
        for key,value in temp.items():
            if key.split("-")[2]=="ref":
                config_dir[f"step{key.split('-')[1]}_ref_species"]=value.split()[0]
                config_dir[f"step{key.split('-')[1]}_ref_release"]=value.split()[1]
            else:
                config_dir[f"step{key.split('-')[1]}_{key.split('-')[2]}"]=value
    config_string = "\n".join(f"{key}: {value}" for key, value in config_dir.items()) + "\n"
    config_string += "samples:\n" + "\n".join(f"  - {sample}" for sample in samples_list) + "\n"
    config_string += "\n".join(f"step{i}_params : \"{params}\"" for i, params in enumerate(params_list)) + "\n"
       
    output_config_list = []
    for i, soft in enumerate(soft_list):
        config = tools_config.get(soft, tools_config[type_list[i]])[mode_list[i]]
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
    config_string += "steps:\n"
    for i in range(len(soft_list)):
        config_string += "  step" + str(i) + ":\n" + "    input:\n"
        config_string += "\n".join(f"      {key} : \"{value}\"" for key, value in input_config_list[i].items()) + "\n"
        config_string += "    output:\n"
        config_string += "\n".join(f"      {key} : \"{value}\"" for key, value in output_config_list[i].items()) + "\n"


 
    smk_dir = paths_dict["smk_dir"]
    snakemake_string = ''.join(f'include: "{os.path.join(smk_dir, filename)}"\n' for filename in os.listdir(smk_dir))
    snakemake_string += "configfile: \"config.yaml\"\n" + "rule all:\n" + "    input:\n"
    snakemake_string += "\n".join(f"        expand(\"{temp}\", sample=config[\"samples\"], output_dir=config[\"output_dir\"])," 
        for temp in [list(record.values())[0] for record in output_config_list]) + "\n"
    return config_string, snakemake_string


def generate_present_script(steps, mode, select_raw_folder, samples_list, task_id, paths_dict):
    input_list = [record["select_option"] for record in steps]  # 这是字典
    soft_list = [record["type"].split()[1] for record in steps]
    types_list = [record["type"] for record in steps]
    params_list = [record["params"] for record in steps]
    mode_list = [record["mode"] for record in steps]

    formatted_string = ""
    for i, type in enumerate(types_list):
        values_dict = {}
        for key, value in input_list[i].items():
            if value.startswith("raw"):
                values_dict[key] = select_raw_folder + "/${sample}" + value.split("raw", 1)[1]
            else:
                prefix = str(value).split("_")[0]
                suffix = value[len(prefix) + 1:].split(".")[0]
                values_dict[key] = f"$output_dir/${{sample}}/{prefix}_{suffix}/{suffix}_${{sample}}" + "" \
                #                   value.split(".raw", 1)[1]

        values_dict["params"] = params_list[i]
        values_dict["step"] = "step" + str(i)
        
        formatted_string += my_code_script.get(type,my_code_script["BAM remove_duplicates"])[mode_list[i]].format(**values_dict)
     

    my_code = (f"mode=\"{mode}\"\n"
               f"output_dir=\"/data/output/{task_id}\"\n")
    my_code += '\n'.join([f'{soft}="{paths_dict.get(soft,soft)}"' for soft in soft_list]) + "\n"
    my_code += 'samples=("{}")'.format('" "'.join(samples_list)) + "\n"
    my_code += 'mkdir -p $output_dir/logs\n'
    my_code += 'for sample in "${samples[@]}"; do\n'
    my_code += formatted_string + "\n"
    my_code += "done\n"

    return my_code






