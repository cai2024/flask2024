import os
import copy
import yaml
from collections import defaultdict
from .utils_common import SoftStep, SoftStepManager, process_all_options


def get_batch_preset(preset_type,mode,soft_config,index_options_dict):
    preset_dict = {
        "WGS": {
            "baseOptions": {"pair": {"fq1": ["raw_1.fq.gz"], "fq2": ["raw_2.fq.gz"]}, "single": {"fq": ["raw.fq.gz"]}},
            "pair": [
                SoftStep(0, "Trim trim_galore", { }, {}, "", "pair", "raw.fq.gz",{}, 
                         {"fq1": "raw_1.fq.gz", "fq2": "raw_2.fq.gz"}, {}),
                SoftStep(1, "map bwa", {}, {}, "", "pair", "raw.fq.gz", {},
                         {"fq1": "step0_trim_galore_raw_1.fq.gz", "fq2": "step0_trim_galore_raw_2.fq.gz"},{}), ],
            "single": [
                SoftStep(0, "Trim trim_galore", {}, {}, "", "single", "raw.fq.gz", {}, 
                         {"fq": "raw.fq.gz"}, {}),
                SoftStep(1, "map bwa", {}, {}, "", "single", "raw.fq.gz", {},
                         {"fq": "step0_trim_galore_raw.fq.gz"}, {}), ]
        },
    }
    return  init_div_steps2( [step.to_json() for step in preset_dict.get(preset_type, preset_dict["WGS"])[mode]],soft_config,index_options_dict) 




def init_div_batch(steps,soft_config,index_options_dict):
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
       
        if not soft_step.rec_option:
            soft_step.update_details(rec_option={
                key: all_options[key][
                    max(enumerate([int(temp.split("_")[0][4:]) + 1 if temp.split("_")[0][4:].isdigit() else 0
                               for temp in all_options[key]]), key=lambda x: x[1])[0]] for key in key_list})
        soft_step.update_details(
            selectableIndex=index_options_dict.get(step["type"].split()[1], index_options_dict["default"]))
        manager.add_step(soft_step)
    return manager.to_json()
