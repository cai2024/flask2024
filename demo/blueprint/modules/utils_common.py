import os
import yaml
from collections import defaultdict
import copy
import pandas as pd
import csv


# 这里all_options的key是数据类型，值是字符串列表

def get_samples_list(input_dir, mode, input_type="fq.gz"):
    input_list = [f for f in os.listdir(input_dir) if
                  f.endswith(input_type) and os.path.isfile(os.path.join(input_dir, f))]
    if input_type == "fq.gz" and mode == "pair":
        input_list = [f.replace('_1.fq.gz', '.fq.gz') for f in input_list if
                      f.endswith('_1.fq.gz') and f.replace('_1.fq.gz', '_2.fq.gz') in input_list]
    return [f[:-len(input_type)-1] for f in input_list]



def find_matching_key_value(input_dict, target_string):
    """
    找到字典中与目标字符串 split(".")[1:] 结果匹配的值并返回。
    
    :param input_dict: 输入字典，键是字符串，值是字符串。
    :param target_string: 输入的目标字符串。
    :return: 返回第一个匹配的字典项的键和值，或 None。
    """
    # 将目标字符串分割
    target_parts = target_string[target_string.rfind('_raw') + len('_raw'):]    
    # 遍历字典中的每个键值对
    # print(target_parts)
    # print(input_dict)
    for key, value in input_dict.items():
        # 将字典值进行分割
        value_parts = value[value.rfind('{sample}') + len('{sample}'):]
        # 如果匹配则返回
        if target_parts == value_parts:
            return value
    
    # 如果没有找到匹配项
    return None


def process_all_options( steps, all_options, soft_config):
    result = [copy.deepcopy(all_options)]
    for step in steps:
        step_id = step["id"]
        step_type, step_soft = step["type"].split()
        config = soft_config.get(step_soft, soft_config[step_type])[step["mode"]]["output"]
        for k, v in config.items():
            all_options[k].extend(f"step{step_id}_" + os.path.basename(f).replace("_{{sample}}", "_raw")
                                  .format(input_soft=step_soft) for f in v.values())
        result.append(copy.deepcopy(all_options))
    return result
#这个函数决定了展示的效果
# 注意到input的值没有任何意义暂时


# 现在的主要问题是，最终输出的时候，怎么处理同一个软件，有多个同类型的不同文件名的文件
# 举例，假设一个函数，生成了两个不同bam，不妨_1和_2,
# 他们作为输出，可以在smk中写下，"{output_dir}/{sample}/{step}_CGS_spikein_modify/CGS_spikein_modify_{sample}.hmc.all.bed"
# 形式上，在{sample}后加.x.类型的方式来处理，目前可以

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

        # config获取对应的每一步输出，这是个字典
        output_config = {k: v.format(input_step=f"step{i}", input_soft=soft) for inner_dict in config.values() for k, v in inner_dict.items()}
        # 这里key是软件内部key， 这里内部就存在bam1这样的非标准 输入

        # 这里key是变量名，value理论也是要变量值也就是字符串，这个config.values()，是一个对一个的字符串
        # 这里保证了output的字符串和变量名，这里的output会显示在前端

       
        output_config_list.append(output_config)
        config_dict["steps"][f"step{i}"]["output"] = output_config

        suffix_map = {"fq1": "_1.fq.gz", "fq2": "_2.fq.gz", "fq": ".fq.gz"}
        # 这里获取input的逻辑，是从列表某一步step字典中，获取k，这里k是选中选项的key和value，k这里应该包含bam1
        # 这里的问题就在这，这里拿变量名去找，虽然在前面的output也是用的变量名，但是二者的变量名又不一致，姑且可以用变量类别对付用
        # 这里又引出一个问题，如果output的变量类别下，有多个同类别的不同变量，怎么办
        # 比如我调用的某一个txt输出下的，hmc/mc，我用变量名去找，前提是我也用hmc，mc来命名，但是实际我可能用txt来命名，这样就对应不上了
        # 这里是个大问题，同类别多输入，多输出
        # 文件命名上是有差异的
        # 比如hmc/mc后醉上不一样啊

        # step i 的input，整体看是个字典，k是选项卡实际变量名，重点是变量值，变量值是字符串，一个变量名对一个字符串
        # 找字符串，本质是找特征
        # 只要字符串能体现，他出生时的变量名，就能在output_config_list【i】找到唯一的对应变量名，而不是用自定义变量名
        # 但是随意起名，又会影响smk最终的生成



        # bam1  step{step_id}_{step_type}.{sample}.bam
        # 去找对应列表i字典中，某个key和value，首先i是由step_id确定的，这个字典，有很多个输出，也就有很多个变量名，假设
        # {sample}.1.bam {sample}.2.bam,这里两个bam，那么他的标识性信息，就在sample的后几位，只要这后面能对应上，需要的字符串就对应上了
        # 出现的第一个.后面的内容本质就是标识了
        # 关键变量值 output_config_list[int(v.split("_")[0][4:])]  对应i，v.split(".")[1:]   对应的该字典的每个value都要做split去找一模一样的


        config_dict["steps"][f"step{i}"]["input"] = {
            k: f"{select_raw_folder}/{{sample}}{suffix_map.get(k, f'.{k}')}"
            if v.startswith("raw")
            else find_matching_key_value(output_config_list[int(v.split("_")[0][4:])], v)  
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
        temp_script = code_script.get(record["type"], code_script["default"])[record["mode"]].replace('{step}', values_dict
                                                                                                          ["step"])
        formatted_string += temp_script.format(**values_dict) + "\n"
    soft_list = [record["type"].split()[1] for record in steps]
    output_dir =config_dict["output_dir"]
    ref_path =config_dict["ref_path"]

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
            self.baseOptions = args[6] if len(args) > 6 else None
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
            self.baseOptions = kwargs.get('baseOptions', None)
            self.selectableOptions = kwargs.get('selectableOptions', None)
            self.rec_option = kwargs.get('rec_option', None)
            self.selectableIndex = kwargs.get('selectableIndex', None)

    def update_details(self, id=None, type=None, select_option=None, select_index=None, params=None, mode=None,
                       baseOptions=None, selectableOptions=None, rec_option=None, selectableIndex=None):
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
        if baseOptions:
            self.baseOptions = baseOptions
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
            'baseOptions': self.baseOptions,
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





def load_config():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(current_dir, 'config.yaml')
    smk_path = os.path.join(current_dir, '../workflows/')
    with open(config_path, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)
        available_steps = config.get('available_steps', {})
        available_steps2 = config.get('available_steps2', {})
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
    return available_steps,available_steps2,soft_config,index_options_dict,paths_dict,code_script

