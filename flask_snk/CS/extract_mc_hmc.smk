rule extract_mc_hmc:
    input:
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/bam/{align_soft}/{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        mc="{output_dir}/bam/{align_soft}/{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=config["ref_genome"],
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """




rule de_extract_mc_hmc:
    input:
        bam="{output_dir}/debam/{align_soft}/de_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/debam/{align_soft}/de_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        mc="{output_dir}/debam/{align_soft}/de_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/de_extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=config["ref_genome"]
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """



rule extract_mc_hmc_sp:
    input:
        bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        fq1="{output_dir}/trim/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir']),
    output:
        hmc= "{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        mc="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
    log:
        "{output_dir}/logs/extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    threads: 20
    params:
        ref=spikein_ref
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type mc ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1}  --ref {params.ref}.fa --parallel {threads} --type hmc ) 2>> {log}
        """

