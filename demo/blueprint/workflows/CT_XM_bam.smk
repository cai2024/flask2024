rule CT_XM_bam:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bam = "{output_dir}/{sample}/{step}_CT_XM_bam/CT_XM_bam_{sample}.bam",
    params:
        ref = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/" +
            config[f"{wildcards.step}_ref_release"] + ".fa"
        ),
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_CT_XM_bam/{sample}.log"
    threads: 1
    shell:
        """
        python3 {params.py_ref}/CT_XM_bam.py --bam {input.bam} --output {output.bam} --ref {params.ref} \
        --parallel {threads}  {params.other} > {log} 2>&1 
       
        """

rule CT_XM_sp_bam:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bam = "{output_dir}/{sample}/{step}_CT_XM_sp_bam/CT_XM_sp_bam_{sample}.bam",
    params:
        ref = config["ref_path"] + "/spikein/spikein.fa",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_CT_XM_sp_bam/{sample}.log"
    threads: 1
    shell:
        """
        python3 {params.py_ref}/CT_XM_bam.py --bam {input.bam} --output {output.bam} --ref {params.ref} \
        --parallel {threads}  {params.other} > {log} 2>&1
        
        """

