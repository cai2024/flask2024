import os
import copy
import yaml
import json
import gzip
import base64
from collections import defaultdict
from io import BytesIO
from graphviz import Digraph
from .utils_common import SoftStep, SoftStepManager, process_all_options






def get_batch_preset(base64_encoded_data,soft_config,index_options_dict):
    compressed_data = base64.urlsafe_b64decode(base64_encoded_data.encode('utf-8'))
    json_data = gzip.decompress(compressed_data).decode('utf-8')
    original_data = json.loads(json_data)
    result=[]
    for temp in original_data:
        result.append(SoftStep(temp[0],temp[1],{},{},temp[2],temp[3],temp[4],{},temp[5],{}))
    return  init_div_batch([step.to_json() for step in result],soft_config,index_options_dict) 




def look_div_batch(steps, soft_config, index_options_dict):
    temp=[]
    for step in steps:
        new_rec_option = {key: step['selectableOptions'][key].index(value) for key, value in step['select_option'].items()}
        temp.append([step["id"],step["type"],step["params"],step["mode"],step["baseOptions"],new_rec_option])
   
    json_data = json.dumps(temp, ensure_ascii=False) 
    compressed_data = gzip.compress(json_data.encode('utf-8'))
    base64_encoded_data = base64.urlsafe_b64encode(compressed_data).decode('utf-8')

    nodes_info = {f"step{k[0]}": k[1] for k in temp}
    nodes_info["raw_data"] = "原始数据"
    data = [["raw_data" if value.startswith("raw") else value.split("_")[0], f"step{k['id']}"] for k in steps for value in k["select_option"].values()]
    graph = Digraph(format='png')
    graph.attr(size='10,10', fontname="Microsoft YaHei")
    for node, label in nodes_info.items():
        graph.node(node, f"{node}\n{label}")
    for source, target in data:
        graph.edge(source, target)
    img_bytes = BytesIO(graph.pipe(format='png'))
    img_base64 = base64.b64encode(img_bytes.getvalue()).decode('utf-8')
   
    return base64_encoded_data,img_base64


def init_div_batch(steps, soft_config, index_options_dict):
    manager = SoftStepManager()
    raw_baseOptions = copy.deepcopy(steps[0]["baseOptions"])
    baseOptions = defaultdict(list)
    baseOptions.update(raw_baseOptions)
    all_options_list = process_all_options(steps, baseOptions, soft_config)
    for i, step in enumerate(steps):
        all_options = all_options_list[i].copy()
        soft_step = SoftStep(**step)
        step_type, step_soft = step["type"].split()
        key_dict = soft_config.get(step_soft, soft_config[step_type])[step["mode"]]["input"]
        soft_step.update_details(selectableOptions={k: all_options.get(key, []).copy() for key, value in key_dict.items() for k in value})
        soft_step.update_details(selectableIndex=index_options_dict.get(step_soft, index_options_dict["default"]))
        if soft_step.params == "":
            soft_step.update_details(params=soft_config.get(step_soft, soft_config[step_type]).get(step["mode"], {}).get("params", ""))
        if not soft_step.rec_option:
            soft_step.update_details(rec_option={
                k: max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                        for temp in all_options[key]]), key=lambda x: x[1])[0] for key, value in key_dict.items() for k in value})
        manager.add_step(soft_step)
    return manager.to_json()

