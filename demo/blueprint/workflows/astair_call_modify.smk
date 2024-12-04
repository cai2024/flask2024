rule astair_call_modify:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bed = "{output_dir}/{sample}/{step}_astair_call_modify/astair_call_modify_{sample}.bed",
    params:
        ref = lambda wildcards: (
            config["ref_path"] + "/" +
            config[f"{wildcards.step}_ref_species"] + "/" +
            config[f"{wildcards.step}_ref_release"] + "/" +
            config[f"{wildcards.step}_ref_release"] + ".fa"
        ),
        out_dir = "{output_dir}/{sample}/{step}_astair_call_modify",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "py39"
    log:
        "{output_dir}/logs/{sample}/{step}_astair_call_modify/{sample}.log"
    threads: 10
    shell:
        """
        astair call -i {input.bam} -f {params.ref}  -t {threads} -d {params.out_dir} {params.other}
        temp_bam=$(basename {input.bam} .bam)_mCtoT_all.mods 
        mv {params.out_dir}/$temp_bam {output.bed}


        """

rule astair_spikein_modify:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        bed = "{output_dir}/{sample}/{step}_astair_spikein_modify/astair_spikein_modify_{sample}.bed",
    params:
        ref = config["ref_path"] + "/spikein/spikein.fa",
        out_dir = "{output_dir}/{sample}/{step}_astair_spikein_modify",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "py39"
    log:
        "{output_dir}/logs/{sample}/{step}_astair_spikein_modify/{sample}.log"
    threads: 10
    shell:
        """
        astair call -i {input.bam} -f {params.ref}  -t {threads} -d {params.out_dir} {params.other}
        temp_bam=$(basename {input.bam} .bam)_mCtoT_all.mods
        mv {params.out_dir}/$temp_bam {output.bed}
        """

