rule CGS_extract_mc_hmc:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"]
    output:
        mc = "{output_dir}/{sample}/{step}_CGS_call_modify/CGS_call_modify_{sample}.mc.all.bed",
        hmc = "{output_dir}/{sample}/{step}_CGS_call_modify/CGS_call_modify_{sample}.hmc.all.bed"
    params:
        ref = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/" +
            config[f"{wildcards.step}_ref_release"] + ".fa"
        ),
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"]
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_CGS_call_modify/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.cut_fq1} --ref {params.ref} \
        --parallel {threads} --type mc --output {output.mc} {params.other} > {log} 2>&1 
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.cut_fq1} --ref {params.ref} \
        --parallel {threads} --type hmc --output {output.hmc} {params.other} > {log} 2>&1
        """
