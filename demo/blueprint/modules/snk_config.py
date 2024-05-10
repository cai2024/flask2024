import os
import subprocess
import yaml

from .utils import get_fq_list
from .utils import ana_get_res


def submit_snakemake(info):
    outdir=info['output_dir']
    if info['ref_ana']=="mouse mm9":
        ref="/data/reference2024/mouse/mm9/bowtie2/mm9"
    elif info['ref_ana']=="mouse mm10":
        ref="/data/reference2024/mouse/mm10/bowtie2/mm10"
    elif info['ref_ana']=="mouse mm39":
        ref="/data/reference2024/mouse/mm30/bowtie2/mm39"
    elif info['ref_ana']=="human hg18":
        ref="/data/reference2024/human/hg18/bowtie2/hg18"
    elif info['ref_ana']=="human hg19":
        ref="/data/reference2024/human/hg19/bowtie2/hg19"
    elif info['ref_ana']=="human hg38":
        ref="/data/reference2024/human/hg38/bowtie2/hg38"
    elif info['ref_ana']=="human hs1":
        ref="/data/reference2024/human/hs1/bowtie2/hs1"

    if info["mode"] == "single":
        fastqs = get_fq_list(info['data_ana'],False)
        sampleList=[temp[:temp.find(".fq.gz")] for temp in fastqs]
    else:
        fastqs = get_fq_list(info['data_ana'])
        sampleList=[temp[:temp.find("_1.fq.gz")] for temp in fastqs]
#---------------------------------------------------------------------------------设置脚本和参数文件----------
    if info['pipeline']=="WGS":
        script_file="/data/cailab/flask2024/flask_snk/WGS.py"
        config = {
            'ref_genome': ref,
            'fq_in_path': info['data_ana'],
            'sampleList': sampleList,
            'fastp_params': info['fastp_params'],
            'output_dir': info['output_dir'] ,
            'mode': info["mode"],
            'bowtie2_params': info['bowtie2_params'],
            'picard_params': info['picard_params']
        }
        with open(outdir+'/config.yaml', 'w') as file:
            yaml.dump(config, file)
    elif info['pipeline']=="CS":
        script_file="/data/cailab/flask2024/flask_snk/CS.py"
        config = {
            'ref_genome': ref,
            'fq_in_path': info['data_ana'],
            'sampleList': sampleList,
            'fastp_params': info['fastp_params'],
            'output_dir': info['output_dir'] ,
            'mode': info["mode"],
            'bowtie2_params': info['bowtie2_params'],
            'picard_params': info['picard_params'],
            'hairpin_cut_params': info['hairpin_cut_params'],
            'spikein_params': info['spikein_params'],
        }
        with open(outdir+'/config.yaml', 'w') as file:
            yaml.dump(config, file)

    else:
        script_file="/data/cailab/flask_snk/11.py"
#------------------------------------------------------------------------------------修改区域-------------------
    snakemake_cmd = [
        'nohup',
        '/data/biosoft/soft2024/conda/anaconda_23.4.7/bin/snakemake',
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



