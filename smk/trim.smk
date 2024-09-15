rule trim_galore_se:
    input:
        fq = lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        trim_fq="{output_dir}/{sample}/{step}_trim_galore/trim_galore_{sample}.fq.gz",
    params:
        trim_dir="{output_dir}/{sample}/{step}_trim_galore",
        other=lambda wildcards: config[f"{wildcards.step}_params"],
        trim_galore=config["trim_galore"],      
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024",  
    log:
        "{output_dir}/logs/{sample}/{step}_trim_galore/trim_galore_{sample}_single.log",
    threads: 3,
    shell:
        """
        {params.trim_galore}  {params.other} --cores {threads} -o {params.trim_dir} {input.fq} > {log} 2>&1
        temp_name=$(basename {input.fq} | sed 's/\.fq\.gz$/_trimmed.fq.gz/')
        mv {params.trim_dir}/$temp_name {output.trim_fq}

                 
        original_count=$(zcat {input.fq} | echo $((`wc -l`/4)))
        trimmed_count=$(zcat {output.trim_fq} | echo $((`wc -l`/4)))
        clean_read_percentage=$(python3 -c "print(format(100 * $trimmed_count / $original_count, '.2f'))")
        echo "Clean_Read_Percentage $clean_read_percentage%" >> {params.final_log}
        """
rule trim_galore_pe:
    input:
        fq1=lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2=lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"],
    output:
        trim_fq1="{output_dir}/{sample}/{step}_trim_galore/trim_galore_{sample}_1.fq.gz",
        trim_fq2="{output_dir}/{sample}/{step}_trim_galore/trim_galore_{sample}_2.fq.gz",
    params:
        trim_dir="{output_dir}/{sample}/{step}_trim_galore",
        other=lambda wildcards: config[f"{wildcards.step}_params"],
        trim_galore=config["trim_galore"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024",  
    log:
        "{output_dir}/logs/{sample}/{step}_trim_galore/trim_galore_{sample}_pair.log",
    threads: 3,
    shell:
        """
        {params.trim_galore}  {params.other} --paired --cores {threads} -o {params.trim_dir} {input.fq1} {input.fq2} > {log} 2>&1
        temp_name1=$(basename {input.fq1} | sed 's/\.fq\.gz$/_val_1.fq.gz/')
        temp_name2=$(basename {input.fq2} | sed 's/\.fq\.gz$/_val_2.fq.gz/')
        mv {params.trim_dir}/$temp_name1 {output.trim_fq1}
        mv {params.trim_dir}/$temp_name2 {output.trim_fq2}

        original_count=$(zcat {input.fq1} | echo $((`wc -l`/4)))
        trimmed_count=$(zcat {output.trim_fq1} | echo $((`wc -l`/4)))
        clean_read_percentage=$(python3 -c "print(format(100 * $trimmed_count / $original_count, '.2f'))")
        echo "Clean_Read_Percentage $clean_read_percentage%" >> {params.final_log}
        """

rule fastp_se:
    input:
        fq=lambda wildcards: config["steps"][wildcards.step]["input"]["fq"]
    output:
        trim_fq="{output_dir}/{sample}/{step}_fastp/fastp_{sample}.fq.gz",
    params:
        json="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_single.json",
        html="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_single.html",
        other=lambda wildcards: config[f"{wildcards.step}_params"],
        fastp=config["fastp"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024",
    log:
        "{output_dir}/logs/{sample}/{step}_fastp/fastp_{sample}_single.log",
    threads: 3,
    shell:
        """
        {params.fastp} {params.other} -i {input.fq} -o {output.trim_fq} -j {params.json} -h {params.html} -w {threads}  > {log} 2>&1

        original_count=$(zcat {input.fq} | echo $((`wc -l`/4)))
        trimmed_count=$(zcat {output.trim_fq} | echo $((`wc -l`/4)))
        clean_read_percentage=$(python3 -c "print(format(100 * $trimmed_count / $original_count, '.2f'))")
        echo "Clean_Read_Percentage $clean_read_percentage%" >> {params.final_log}
        """
rule fastp_pe:
    input:
        fq1=lambda wildcards: config["steps"][wildcards.step]["input"]["fq1"],
        fq2=lambda wildcards: config["steps"][wildcards.step]["input"]["fq2"],
    output:
        trim_fq1="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_1.fq.gz",
        trim_fq2="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_2.fq.gz",
    params:
        json="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_pair.json",
        html="{output_dir}/{sample}/{step}_fastp/fastp_{sample}_pair.html",
        other=lambda wildcards: config[f"{wildcards.step}_params"],
        fastp=config["fastp"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024",
    log:
        "{output_dir}/logs/{sample}/{step}_fastp/fastp_{sample}_pair.log",
    threads: 3,
    shell:
        """
        {params.fastp} {params.other} -i {input.fq1} -I {input.fq2} -o {output.trim_fq1} -O {output.trim_fq2} -j {params.json} -h {params.html} -w {threads} > {log} 2>&1

        original_count=$(zcat {input.fq1} | echo $((`wc -l`/4)))
        trimmed_count=$(zcat {output.trim_fq1} | echo $((`wc -l`/4)))
        clean_read_percentage=$(python3 -c "print(format(100 * $trimmed_count / $original_count, '.2f'))")
        echo "Clean_Read_Percentage $clean_read_percentage%" >> {params.final_log}
        """


