
if config["mode"] == "single":
    se = "-se"
else:
    se = ""


rule process_bam:
    input:
        bam="{output_dir}/{{sample}}/de_bam/{dup_soft}/{{sample}}_dedup.bam".format(output_dir=config['output_dir'],dup_soft=config['dup_soft'])
    output:
        bam="{output_dir}/{{sample}}/de_bam/merge/de_{{sample}}.bam".format(output_dir=config['output_dir'])
    conda:
        "flask2024",
    log:
        "{output_dir}/{{sample}}/logs/process_genome/de_{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    shell:
        """
        ({bamUtil} clipOverlap --in {input.bam} --out {output.bam} --stats ) 2> {log}
        samtools index {output.bam}
        """


rule astair_call_CpG: 
    input:
        bam="{output_dir}/{{sample}}/de_bam/merge/de_{{sample}}.bam".format(output_dir=config['output_dir'])
    output:
        "{output_dir}/{{sample}}/de_bam/{call_soft}/de_{{sample}}_mCtoT_CpG.mods.gz".format(output_dir=config['output_dir'],call_soft=config['call_soft'])
    conda:
        "py39",
    log:
        "{output_dir}/{{sample}}/logs/call/de_{{sample}}.log".format(output_dir=config['output_dir']),
    params:
        other=config['call_params'],
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] + "/" + config["ref_release"] + ".fa",
        outdir="{output_dir}/{{sample}}/de_bam/{call_soft}".format(output_dir=config['output_dir'],call_soft=config['call_soft'])
    threads: 10
    shell:
        """
        (astair call -sc True  -i {input.bam} -f {params.ref} -t {threads} -m mCtoT --context CpG --gz -d {params.outdir} {params.other} {se}) 2> {log}
        """


rule astair_call_all: 
    input:
        bam="{output_dir}/{{sample}}/de_bam/merge/de_{{sample}}.bam".format(output_dir=config['output_dir'])
    output:
        "{output_dir}/{{sample}}/de_bam/{call_soft}/de_{{sample}}_mCtoT_all.mods.gz".format(output_dir=config['output_dir'],call_soft=config['call_soft'])
    log:
        "{output_dir}/{{sample}}/logs/call/de_{{sample}}_all.log".format(output_dir=config['output_dir']),
    params:
        other=config['call_params'],
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] + "/" + config["ref_release"] + ".fa",
        outdir="{output_dir}/{{sample}}/de_bam/{call_soft}".format(output_dir=config['output_dir'],call_soft=config['call_soft'])
    threads: 10
    conda:
        'py39',
    shell:
        """
        (astair call -sc True -i {input.bam} -f {params.ref} -t {threads} -m mCtoT --context all --gz -d {params.outdir} {params.other} {se} ) 2> {log}
        """




