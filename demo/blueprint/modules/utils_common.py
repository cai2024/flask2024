import os
import yaml
from collections import defaultdict
import copy
import pandas as pd
import csv




def get_fq_list(indir, paired_end=True):
    fastqs = []
    if paired_end:
        for file in os.listdir(indir):
            if file.endswith('1.fq.gz') and file[:-7] + '2.fq.gz' in os.listdir(indir):
                fastqs.append(file[:-7] + "1.fq.gz" + " " + file[:-7] + "2.fq.gz")
            if file.endswith('1.fastq.gz') and file[:-10] + '2.fastq.gz' in os.listdir(indir):
                fastqs.append(file[:-10] + "1.fastq.gz" + " " + file[:-10] + "2.fastq.gz")
    else:
        fastqs = [f for f in os.listdir(indir) if f.endswith('.fq.gz') or f.endswith('.fastq.gz')]
    return fastqs


def get_samples_list(select_raw_folder, fq_type, raw_myOptions):
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




def generate_steps_config(steps, select_raw_folder, samples_list, task_id, soft_config, paths_dict, code_script):
    config_dict = paths_dict.copy()
    config_dict["output_dir"] = f"/data/cailab/flask_out/{task_id}"
    config_dict["samples"] = samples_list
    config_dict.update({key: value for record in steps for key, value in record["select_index"].items()})
    config_dict.update({f"step{i}_params": record["params"] for i, record in enumerate(steps)})
    output_config_list = []
    config_dict["steps"] = defaultdict(lambda: defaultdict(dict))
    for i, step in enumerate(steps):
        type_name, soft = step["type"].split()[0], step["type"].split()[1]
        config = soft_config.get(soft, soft_config.get(type_name))[step["mode"]]["output"]
        output_config = {k: v.format(input_step=f"step{i}", input_soft=soft) for k, v in config.items()}
        output_config_list.append(output_config)
        config_dict["steps"][f"step{i}"]["output"] = output_config
        suffix_map = {"fq1": "_1.fq.gz", "fq2": "_2.fq.gz", "fq": ".fq.gz"}
        config_dict["steps"][f"step{i}"]["input"] = {
            k: f"{select_raw_folder}/{{sample}}{suffix_map.get(k, f'.{k}')}"
            if v.startswith("raw")
            else output_config_list[int(v.split("_")[0][4:])][k]
            for k, v in step["select_option"].items()
        }
    config_dict["steps"] = {k: dict(v) for k, v in config_dict["steps"].items()}
    config_string = yaml.dump(dict(config_dict), default_flow_style=False)
    smk_dir = paths_dict["smk_dir"]
    snakemake_string = ''.join(f'include: "{os.path.join(smk_dir, filename)}"\n' for filename in os.listdir(smk_dir) if filename.endswith('.smk') )
    snakemake_string += "configfile: \"config.yaml\"\n" + "rule all:\n" + "  input:\n"
    snakemake_string += "\n".join(('    expand("{output}", sample=config["samples"], output_dir=config["output_dir"]),'
                                   .format(output=list(temp.values())[0])) for temp in output_config_list) + "\n"

    formatted_string = ""
    for i, record in enumerate(steps):
        values_dict = {}
        for key, value in record["select_option"].items():
            if value.startswith("raw"):
                values_dict[key] = select_raw_folder + "/${sample}" + value.split("raw", 1)[1]
            else:
                prefix = str(value).split("_")[0]
                suffix = value[len(prefix) + 1:].split(".")[0]
                values_dict[key] = f"$output_dir/${{sample}}/{prefix}_{suffix}/{suffix}_${{sample}}"
        values_dict["params"] = record["params"]
        values_dict["step"] = f"step{i}"
        values_dict.update(record["select_index"])
        temp_script = code_script.get(record["type"], code_script["default"])[record["mode"]].replace('{step}', values_dict["step"])
        formatted_string += temp_script.format(**values_dict) + "\n"
    soft_list = [record["type"].split()[1] for record in steps]
    output_dir=config_dict["output_dir"]    
    ref_path=config_dict["ref_path"]
    
    my_code = '\n'.join([f'{soft}="{paths_dict.get(soft, soft)}"' for soft in soft_list]) + "\n"
    my_code += f"ref_path=\"{ref_path}\"\n"
    my_code += f"output_dir=\"{output_dir}\"\n"
    my_code += 'samples=("{}")'.format('" "'.join(samples_list)) + "\n"
    my_code += 'mkdir -p $output_dir/logs\n'
    my_code += 'for sample in "${samples[@]}"; do\n' + formatted_string + "\n" + "done\n"
    return config_string, snakemake_string, my_code






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

    def update_details(self, id=None, type=None, select_option=None, select_index=None, params=None, mode=None,
                       raw_type=None, selectableOptions=None, rec_option=None, selectableIndex=None):
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


def process_all_options(soft_config, steps, all_options):
    result = [copy.deepcopy(all_options)]
    for step in steps:
        config = soft_config.get(step["type"].split()[1], soft_config[step["type"].split()[0]])[step["mode"]]["output"]
        for key in config.keys():
            temp_string = config[key].split("/")[-1].replace("_{{sample}}", "_raw")
            step_id = step["id"]
            input_soft = step["type"].split()[1]
            formatted_string = f"step{step_id}_" + temp_string.format(input_soft=input_soft)
            all_options[key.split("_")[0]].append(formatted_string)
        result.append(copy.deepcopy(all_options))
    return result



def load_config():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(current_dir, 'config.yaml')
    smk_path = os.path.join(current_dir, '../workflows/')
    with open(config_path, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)
        available_steps = config.get('available_steps', {})
        soft_config = config.get('soft_config', {})
        index_options_dict = config.get('index_options_dict', {})
        paths_dict = config.get('paths_dict', {})
    code_script = {}
    for filename in os.listdir(smk_path):
        if filename.endswith('.yaml'):
            filepath = os.path.join(smk_path, filename)
            with open(filepath, 'r') as file:
                data = yaml.safe_load(file)
                code_script.update(data)
    return available_steps,soft_config,index_options_dict,paths_dict,code_script

