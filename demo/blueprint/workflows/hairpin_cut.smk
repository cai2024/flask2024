
rule hairpin_cut:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        cut_fq = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}.fq.gz",
        cut_fq1 = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}_1.fq.gz",
        cut_fq2 = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}_2.fq.gz"
    params:        
        cut_fq = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}.fq",
        cut_fq1 = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}_1.fq",
        cut_fq2 = "{output_dir}/{sample}/{step}_hairpin_cut/hairpin_cut_{sample}_2.fq",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_hairpin_cut/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/new_hairpin_cut.py --fq1 {input.fq1} --fq2 {input.fq2} --cut_fq1 {params.cut_fq1} \
        --cut_fq2 {params.cut_fq2} --cut_fq {params.cut_fq} {params.other} --parallel {threads} > {log} 2>&1 
        python3 {params.py_ref}/wild_stats.py --fq1 {params.cut_fq1} --fq2 {params.cut_fq2} >> {log} 2>&1
        gzip {params.cut_fq1} {params.cut_fq2} {params.cut_fq}

        temp=$(grep 'resolved ratio:' {log} | sed 's/resolved ratio: /match_model_ratio /g' | head -n 1 )
        echo "{wildcards.step} $temp" >> {params.final_log}
        """


rule hairpin_cut_N:
    input:
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"]
    output:
        cut_fq = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}.fq.gz",
        cut_fq1 = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}_1.fq.gz",
        cut_fq2 = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}_2.fq.gz"
    params:
        cut_fq = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}.fq",
        cut_fq1 = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}_1.fq",
        cut_fq2 = "{output_dir}/{sample}/{step}_hairpin_cut_N/hairpin_cut_N_{sample}_2.fq",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_hairpin_cut/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/new_hairpin_cut_N.py --fq1 {input.fq1} --fq2 {input.fq2} --cut_fq1 {params.cut_fq1} \
        --cut_fq2 {params.cut_fq2} --cut_fq {params.cut_fq} {params.other} --parallel {threads} > {log} 2>&1
        python3 {params.py_ref}/wild_stats.py --fq1 {params.cut_fq1} --fq2 {params.cut_fq2} >> {log} 2>&1
        gzip {params.cut_fq1} {params.cut_fq2} {params.cut_fq}

        temp=$(grep 'resolved ratio:' {log} | sed 's/resolved ratio: /match_model_ratio /g' | head -n 1 )
        echo "{wildcards.step} $temp" >> {params.final_log}
        """

