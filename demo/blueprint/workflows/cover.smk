
rule qualimap_bamqc:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"]
    output:
        cover = "{output_dir}/{sample}/{step}_qualimap_cover/qualimap_cover_{sample}.txt"
    params:
        out_dir = "{output_dir}/{sample}/{step}_qualimap_cover",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        qualimap = config["qualimap"],
        py_ref = config["py_ref"],
        final_log="{output_dir}/logs/{sample}.log",
    conda:
        "flask2024"
    log:
        "{output_dir}/logs/{sample}/{step}_qualimap_cover/{sample}.log"
    threads: 10
    shell:
        """
        {params.qualimap} bamqc -bam {input.bam} -outdir {params.out_dir} -outformat PDF:HTML -nt {threads} --java-mem-size=10G {params.other} > {log} 2>&1 
        cp {params.out_dir}/genome_results.txt {output.cover}


        cover_depth=$(python3 {params.py_ref}/get_log.py cover_depth --input_file {output.cover} )
        var1=$(python3 -c "print('$cover_depth'.split('\t')[0])")
        var2=$(python3 -c "print('$cover_depth'.split('\t')[1])")
        echo "{wildcards.step} cover_ratio $var1" >> {params.final_log}
        echo "{wildcards.step} average_depth $var2" >> {params.final_log}
        """
