import os
import copy
import yaml
from collections import defaultdict
from .utils_common import SoftStep, SoftStepManager, process_all_options



def batch_preset(presetType, mode, soft_config,index_options_dict,batch_preset_dict):
    manager = SoftStepManager()
    preset_list = batch_preset_dict[presetType][mode]
    default_dict = defaultdict(list, batch_preset_dict[presetType]["baseOptions"][mode])
    all_options_list = process_all_options(soft_config, [soft_step.to_json() for soft_step in preset_list],
                                           default_dict)
    for i, soft_step in enumerate(preset_list):
        soft_step.update_details(selectableOptions=all_options_list[i])
        soft_step.update_details(
            selectableIndex=index_options_dict.get(soft_step.type.split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()


def init_div_steps2(steps,soft_config,index_options_dict):
    manager = SoftStepManager()
    baseOptions = defaultdict(list)

    if steps[0]["mode"] == "pair" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq1"] = ["raw_1.fq.gz"]
        baseOptions["fq2"] = ["raw_2.fq.gz"]
    elif steps[0]["mode"] == "single" and steps[0]["raw_type"] == "raw.fq.gz":
        baseOptions["fq"] = ["raw.fq.gz"]
    else:
        baseOptions["bam"] = ["raw.bam"]

    all_options_list = process_all_options(soft_config, steps, baseOptions)

    for i, step in enumerate(steps):
        all_options = all_options_list[i].copy()
        soft_step = SoftStep(**step)
        key_list = soft_config.get(step["type"].split()[1], soft_config[step["type"].split()[0]])[step["mode"]]["input"].keys()

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
