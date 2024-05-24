import os
import subprocess
import yaml

from .utils import get_fq_list
from .utils import ana_get_res


def submit_snakemake(info):
    script_path="/data/cailab/flask2024/flask_snk/"
    snakemake="/data/biosoft/soft2024/conda/anaconda_23.4.7/bin/snakemake"
    if info["align_soft"]=="bowtie2":
        ref="/data/reference2024/" + info['ref_ana'].replace(" ","/")  +"/bowtie2/" + info['ref_ana'].split(" ")[1]
    elif info["align_soft"]=="bismark":
        ref="/data/reference2024/" + info['ref_ana'].replace(" ","/")  +"/bismark/" 



    config = info.copy()
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



