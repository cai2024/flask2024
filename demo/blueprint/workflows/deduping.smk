rule remove_duplicates:
    input:
        bam = lambda wildcards: config["steps"][wildcards.step]["input"]["bam"],
    output:
        de_bam = "{output_dir}/{sample}/{step}_remove_duplicates/remove_duplicates_{sample}.bam"
    log:
        "{output_dir}/logs/{sample}/{step}_remove_duplicates/{sample}.log"
    threads: 4
    conda:
        "flask2024",
    params:
        tmp_report="{output_dir}/{sample}/{step}_remove_duplicates/{sample}.txt",
        other = lambda wildcards: config[f"{wildcards.step}_params"],
        java = config["java"],
        picard = config["picard"],
        final_log="{output_dir}/logs/{sample}.log",
    shell:
        """
        ({params.java} -jar {params.picard} MarkDuplicates --REMOVE_DUPLICATES {params.other} -I {input.bam} -O {output.de_bam} -M {params.tmp_report}) > {log} 2>&1
        samtools index {output.de_bam}    
   
        dedup_reads=$(samtools view -c -F 4 {output.de_bam})
        total_reads=$(samtools view -c -F 4 {input.bam})
        duplication_percentage=$(python3 -c "print(format(100 - (100 * $dedup_reads / $total_reads), '.2f'))" )      
        echo "{wildcards.step} duplication_percentage $duplication_percentage%" >> {params.final_log}
        """

