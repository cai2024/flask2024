rule extract_mc_hmc:
    input:
        bam="{output_dir}/bam/{align_soft}/{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        cut_fq1="{output_dir}/trim/{cut_soft}/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        hmc= "{output_dir}/bam/CGS_call_modify/{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/bam/CGS_call_modify/{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    conda:
        "flask2024",
    threads: 20
    params:
        ref=config["ref_fa"],
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type mc --output {output.mc} ) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type hmc --output {output.hmc}) 2>> {log}
        """


rule de_extract_mc_hmc:
    input:
        de_bam="{output_dir}/de_bam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft']),
        cut_fq1="{output_dir}/trim/{cut_soft}/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        hmc= "{output_dir}/de_bam/CGS_call_modify/{{sample}}_dedup.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/de_bam/CGS_call_modify/{{sample}}_dedup.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/de_extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    conda:
        "flask2024",
    threads: 20
    params:
        ref=config["ref_fa"]
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.de_bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type mc --output {output.mc}) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.de_bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type hmc --output {output.hmc}) 2>> {log}
        """



rule extract_mc_hmc_sp:
    input:
        sp_bam="{output_dir}/sp_bam/{align_soft}/sp_{{sample}}.bam".format(output_dir=config['output_dir'],align_soft=config['align_soft']),
        cut_fq1="{output_dir}/trim/{cut_soft}/{{sample}}_cut_f1.fq".format(output_dir=config['output_dir'],cut_soft=config['cut_soft']),
    output:
        hmc= "{output_dir}/sp_bam/CGS_call_modify/sp_{{sample}}.bam.all.bed.hmc".format(output_dir=config['output_dir']),
        mc="{output_dir}/sp_bam/CGS_call_modify/sp_{{sample}}.bam.all.bed.mc".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/logs/sp_extract_mc_hmc/{{sample}}.log".format(output_dir=config['output_dir'])
    conda:
        "flask2024",
    threads: 20
    params:
        ref=f"{ref_path}/spikein/spikein.fa",
    shell:
        """
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.sp_bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type mc --output {output.mc}) 2> {log}
        (python3 {py_ref}/CGS_extract_modify.py --bam {input.sp_bam} --cutfq1 {input.cut_fq1}  --ref {params.ref} --parallel {threads} --type hmc --output {output.hmc}) 2>> {log}
        """

