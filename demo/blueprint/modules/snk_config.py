import os
import subprocess
import yaml
from .paths import paths_dict
from .utils import get_fq_list
from .utils import ana_get_res

def set_reference_paths(info):
    ref_path=paths_dict["ref_path"]
    ref = ""
    ref_fa = ""
    if info["align_soft"] == "bowtie2":
        ref = ref_path+"/" + info['ref_ana'].replace(" ", "/") + "/bowtie2/" + info['ref_ana'].split(" ")[1]
        ref_fa = ref_path + "/" + info['ref_ana'].replace(" ", "/") + "/" + info['ref_ana'].split(" ")[1] + ".fa"
    elif info["align_soft"] == "bismark":
        ref = ref_path+"/" + info['ref_ana'].replace(" ", "/") + "/bismark/"
        ref_fa = ref + info['ref_ana'].split(" ")[1] + ".fa"
    elif info["align_soft"] == "bwa":
        ref = ref_path + "/"  + info['ref_ana'].replace(" ", "/") + "/bwa/" + info['ref_ana'].split(" ")[1] + ".fa"
        ref_fa = ref    
    return ref, ref_fa


def submit_snakemake(info):
    config = info.copy()
    config.update(paths_dict)
    script_path=config["py_ref"][:-3]
    snakemake=config["snakemake"]
    ref, ref_fa = set_reference_paths(info)

    config['ref_fa']=ref_fa
    config['ref_genome']=ref
    fastqs = get_fq_list(info['fq_in_path'], info["mode"] != "single")
    suffix = ".fq.gz" if info["mode"] == "single" else "_1.fq.gz"
    sampleList = [temp[:temp.find(suffix)] for temp in fastqs]
    config['sampleList']=sampleList



    script_file= script_path + info['pipeline'] + "/"+ info['pipeline'] + ".py"
    with open(info['output_dir']+'/config.yaml', 'w') as file:
        yaml.dump(config, file)




#------------------------------------------------------------------------------------修改区域-------------------
    outdir=info['output_dir']
    snakemake_cmd = [
        'nohup',
        
        snakemake,
        '--snakefile', script_file,
        '-c', str(info['cpu_count']),
        '--use-conda',
        '--conda-frontend', 'conda',
        '--configfile', outdir+'/config.yaml',
        '--default-resources',  f"tmpdir='{outdir}/temp_test'",
    ]

    try:
        with open(outdir + '/snk.log', 'w') as f:
            subprocess.run(snakemake_cmd, check=True, stdout=f, stderr=subprocess.STDOUT)
        ana_get_res(outdir+"/final_log")
    except subprocess.CalledProcessError as e:
        with open(outdir + '/snk.log', 'a') as f:
            print(f"运行出错:{e}",file=f)
            print("请确认参数是否正确，或者@wangzc", file = f)
    return None



