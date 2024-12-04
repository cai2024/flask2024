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
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_CGS_call_modify/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1} --ref {params.ref} \
        --parallel {threads} --type mc --output {output.mc} {params.other} > {log} 2>&1 
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1} --ref {params.ref} \
        --parallel {threads} --type hmc --output {output.hmc} {params.other} > {log} 2>&1

        hmc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.hmc}) 
        hmc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.hmc}) 
        hmc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.hmc})
        mc_cpg=$(awk '$7 == "CpG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.mc})
        mc_chg=$(awk '$7 == "CHG" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.mc}) 
        mc_chh=$(awk '$7 == "CHH" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {output.mc})
  
        hmc=$(python3 -c "print(format($hmc_cpg, '.4f') + '|' + format($hmc_chg, '.4f') + '|' + format($hmc_chh, '.4f'))") 
        mc=$(python3 -c "print(format($mc_cpg, '.4f') + '|' + format($mc_chg, '.4f') + '|' + format($mc_chh, '.4f'))")
        echo "{wildcards.step} hmc_level(CpG|CHG|CHH)  $hmc" >> {params.final_log}
        echo "{wildcards.step} mc_level(CpG|CHG|CHH)  $mc" >> {params.final_log}
        """

rule CGS_extract_sp_mc_hmc:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
        fq1 = lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"]
    output:
        mc = "{output_dir}/{sample}/{step}_CGS_spikein_modify/CGS_spikein_modify_{sample}.mc.all.bed",
        hmc = "{output_dir}/{sample}/{step}_CGS_spikein_modify/CGS_spikein_modify_{sample}.hmc.all.bed",
    params:
        ref = config["ref_path"] + "/spikein/spikein.fa",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        py_ref = config["py_ref"],
        bedtools = config["bedtools"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_CGS_spikein_modify/{sample}.log"
    threads: 10
    shell:
        """
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1} --ref {params.ref} \
        --parallel {threads} --type mc --output {output.mc} {params.other} > {log} 2>&1
        python3 {params.py_ref}/CGS_extract_modify.py --bam {input.bam} --cutfq1 {input.fq1} --ref {params.ref} \
        --parallel {threads} --type hmc --output {output.hmc} {params.other} > {log} 2>&1
       
 
        all_spikein=$(python3 {params.py_ref}/get_log.py spikein --bam {input.bam} --hmc {output.hmc} --mc {output.mc} --bedtools {params.bedtools})
        echo "{wildcards.step} lambda_(cover|dep) $(echo "$all_spikein" | cut -f1)" >> {params.final_log}
        echo "{wildcards.step} lambda_hmc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f2)" >> {params.final_log}
        echo "{wildcards.step} lambda_mc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f3)" >> {params.final_log}
        echo "{wildcards.step} mc_(cover|dep) $(echo "$all_spikein" | cut -f4)" >> {params.final_log}
        echo "{wildcards.step} mc_hmc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f5)" >> {params.final_log}
        echo "{wildcards.step} mc_mc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f6)" >> {params.final_log}
        echo "{wildcards.step} G5hmc_(cover|dep) $(echo "$all_spikein" | cut -f7)" >> {params.final_log}
        echo "{wildcards.step} G5hmc_hmc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f8)" >> {params.final_log}
        echo "{wildcards.step} G5hmc_mc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f9)" >> {params.final_log}
        echo "{wildcards.step} puc19_(cover|dep) $(echo "$all_spikein" | cut -f10)" >> {params.final_log}
        echo "{wildcards.step} puc19_hmc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f11)" >> {params.final_log}
        echo "{wildcards.step} puc19_mc(CpG|CHG|CHH) $(echo "$all_spikein" | cut -f12)" >> {params.final_log}
        """

