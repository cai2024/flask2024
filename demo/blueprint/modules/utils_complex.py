import os
import copy
from collections import defaultdict

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
    "input":{
        "pair":{},
        "single":{}
    },

}
tools_config = {
    "input":{
        "pair":{},
        "single":{}
    },
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
        "ref": {"mouse": ["mm9", "mm10", "mm39"],
                "human": ["hg18", "hg19", "hg38", "hs1"]},
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



def generate_steps_config(steps, select_raw_folder, samples_list, task_id, paths_dict):
    index_list = [record["select_index"] for record in steps if record["select_index"]]
    input_list = [record["select_option"] for record in steps]
    soft_list = [record["type"].split()[1] for record in steps]
    type_list = [record["type"].split()[0] for record in steps]
    params_list = [record["params"] for record in steps]
    mode_list = [record["mode"] for record in steps]
    config_dir = paths_dict.copy()
    config_dir["output_dir"] = f"/data/cailab/flask_out/{task_id}"
    config_dir.update({ key: value for i, temp in enumerate(index_list) for key, value in temp.items() })
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
    snakemake_string += "\n".join(
        f"        expand(\"{temp}\", sample=config[\"samples\"], output_dir=config[\"output_dir\"]),"
        for temp in [list(record.values())[0] for record in output_config_list]) + "\n"
    return config_string, snakemake_string




def generate_present_script(steps, select_raw_folder, samples_list, task_id, paths_dict):
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
                values_dict[key] = f"$output_dir/${{sample}}/{prefix}_{suffix}/{suffix}_${{sample}}" + ""
        values_dict["params"] = params_list[i]
        values_dict["step"] = "step" + str(i)
        formatted_string += my_code_script.get(type, my_code_script["BAM remove_duplicates"])[mode_list[i]].format(**values_dict)

    my_code = f"output_dir=\"/data/output/{task_id}\"\n"
    my_code += '\n'.join([f'{soft}="{paths_dict.get(soft, soft)}"' for soft in soft_list]) + "\n"
    my_code += 'samples=("{}")'.format('" "'.join(samples_list)) + "\n"
    my_code += 'mkdir -p $output_dir/logs\n'
    my_code += 'for sample in "${samples[@]}"; do\n' + formatted_string + "\n" + "done\n"
    return my_code


class SoftStep:
    def __init__(self, *args, **kwargs):
        if args:
            self.id = args[0] if len(args) > 0 else None
            self.type = args[1] if len(args) > 1 else None
            self.select_option = args[2] if len(args) > 2 else None
            self.select_index = args[3] if len(args) > 3 else None
            self.params = args[4] if len(args) > 4 else None
            self.mode = args[5] if len(args) > 5 else None
            self.raw_type = args[6] if len(args) > 6 else None
            self.selectableOptions = args[7] if len(args) > 7 else None
            self.rec_option = args[8] if len(args) > 8 else None
            self.selectableIndex = args[9] if len(args) > 9 else None
        else:
            self.id = kwargs.get('id', None)
            self.type = kwargs.get('type', None)
            self.select_option = kwargs.get('select_option', None)
            self.select_index = kwargs.get('select_index', None)
            self.params = kwargs.get('params', None)
            self.mode = kwargs.get('mode', None)
            self.raw_type = kwargs.get('raw_type', None)
            self.selectableOptions = kwargs.get('selectableOptions', None)
            self.rec_option = kwargs.get('rec_option', None)
            self.selectableIndex = kwargs.get('selectableIndex', None)
    def update_details(self, id=None, type=None, select_option=None,select_index=None, params=None, mode=None, raw_type=None, selectableOptions=None,rec_option=None, selectableIndex=None):
        if id:
            self.id = id
        if type:
            self.type = type
        if select_option:
            self.select_option = select_option
        if select_index:
            self.select_index = select_index
        if params:
            self.params = params
        if mode:
            self.mode = mode
        if raw_type:
            self.raw_type = raw_type
        if selectableOptions:
             self.selectableOptions = selectableOptions
        if rec_option:
            self.rec_option = rec_option
        if selectableIndex:
             self.selectableIndex = selectableIndex

    def to_json(self):
        return {
            'id': self.id,
            'type': self.type,
            'select_option': self.select_option,
            'select_index': self.select_index,
            'params': self.params,
            'mode': self.mode,
            'raw_type': self.raw_type,
            'selectableOptions': self.selectableOptions,
            'rec_option': self.rec_option,
            'selectableIndex': self.selectableIndex,
        }


class SoftStepManager:
    def __init__(self):
        self.steps = []

    def add_step(self, step):
        self.steps.append(step)

    def to_json(self):
        return [step.to_json() for step in self.steps]

batch_preset_dict = {
    "WGS": {
        "baseOptions": {"pair":{"fq1":["raw_1.fq.gz"],"fq2":["raw_2.fq.gz"]},"single":{"fq":["raw.fq.gz"]}},
        "pair": [SoftStep(0,"Trim trim_galore",{"fq1":"raw_1.fq.gz","fq2":"raw_2.fq.gz"},{},"","pair","raw.fq.gz", {},{"fq1":"raw_1.fq.gz","fq2":"raw_2.fq.gz"},{}),
                 SoftStep(1,"map bwa",{"fq1":"step0_trim_galore_raw_1.fq.gz","fq2":"step0_trim_galore_raw_2.fq.gz"},{"ref":"mouse mm9"},"","pair","raw.fq.gz", {},{"fq1":"step0_trim_galore_raw_1.fq.gz","fq2":"step0_trim_galore_raw_2.fq.gz"},{}),],
        "single": [SoftStep(0,"Trim trim_galore",{"fq":"raw.fq.gz"},{},"","single","raw.fq.gz", {},{"fq":"raw.fq.gz"},{}),
                   SoftStep(1,"map bwa",{"fq":"step0_trim_galore_raw.fq.gz"},{"ref":"mouse mm9"},"","single","raw.fq.gz", {},{"fq":"step0_trim_galore_raw.fq.gz"},{}),]
    },
}

def batch_preset(presetType, mode):
    manager = SoftStepManager()
    preset_list = batch_preset_dict[presetType][mode]
    default_dict = defaultdict(list, batch_preset_dict[presetType]["baseOptions"][mode])
    all_options_list = process_all_options(tools_config, [soft_step.to_json() for soft_step in preset_list], default_dict)
    for i,soft_step in enumerate(preset_list):
        soft_step.update_details(selectableOptions=all_options_list[i])
        soft_step.update_details(selectableIndex=index_options_dict.get(soft_step.type.split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()



def process_all_options(tools_config, steps, all_options):
    result = [copy.deepcopy(all_options)]
    for step in steps:
        config = tools_config.get(step["type"].split()[1], tools_config[step["type"].split()[0]])[step["mode"]]
        for key in config.keys():
            temp_string = config[key].split("/")[-1].replace("_{{sample}}", "_raw")
            step_id = step["id"]
            input_soft = step["type"].split()[1]
            formatted_string = f"step{step_id}_" + temp_string.format(input_soft=input_soft)
            all_options[key.split("_")[0]].append(formatted_string)
        result.append(copy.deepcopy(all_options))
    return  result

def init_div_steps2(steps):
    manager = SoftStepManager()
    baseOptions = defaultdict(list)

    if steps[0]["mode"] == "pair" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq1"] = ["raw_1.fq.gz"]
        baseOptions["fq2"] = ["raw_2.fq.gz"]
    elif steps[0]["mode"] == "single" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq"] = ["raw.fq.gz"]
    else:
        baseOptions["bam"] = ["raw.bam"]

    all_options_list = process_all_options(tools_config, steps, baseOptions)

    for i,step in enumerate(steps):
        all_options = all_options_list[i].copy()
        soft_step = SoftStep(**step)
        key_list = tools_input.get(step["type"].split()[1], tools_input[step["type"].split()[0]])[step["mode"]].keys()

        soft_step.update_details(
            selectableOptions={key: all_options.get(key.split("_")[0], []).copy() for key in key_list})
        soft_step.update_details(rec_option={
            key: all_options[key][
                max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                               for temp in all_options[key]]), key=lambda x: x[1])[0]] for key in key_list})
        soft_step.update_details(
            selectableIndex=index_options_dict.get(step["type"].split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()


def init_div_steps3(steps):
    manager = SoftStepManager()
    baseOptions = defaultdict(list, steps[0]["selectableOptions"])
   


    all_options_list = process_all_options(tools_config, steps, baseOptions)

    for i,step in enumerate(steps):
        all_options = all_options_list[i].copy()
        soft_step = SoftStep(**step)
        key_list = tools_input.get(step["type"].split()[1], tools_input[step["type"].split()[0]])[step["mode"]].keys()

        soft_step.update_details(
            selectableOptions={key: all_options.get(key.split("_")[0], []).copy() for key in key_list})
        soft_step.update_details(rec_option={
            key: all_options[key][
                max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                               for temp in all_options[key]]), key=lambda x: x[1])[0]] for key in key_list})
        soft_step.update_details(
            selectableIndex=index_options_dict.get(step["type"].split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()

