rule star_se: 
    input:
        trim_fq1="{output_dir}/{{sample}}/trim/{trim_soft}/{{sample}}_trimmed.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft'])
    output:
        bam="{output_dir}/{{sample}}/bam/star/{{sample}}_single.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/{{sample}}/logs/star/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] +"/star/"+config["ref_release"]+"_"+config["ref_len"],
        bam="{output_dir}/{{sample}}/bam/{{sample}}.star.".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({STAR} --genomeDir {params.ref} --readFilesIn {input.trim_fq1}  --outFileNamePrefix {params.bam} --runThreadN {threads} \
             --outSAMtype BAM Unsorted      --outSAMmultNmax 1 --readFilesCommand zcat {params.other} ) > {log} 2>&1
        
     
        samtools sort {params.bam}Aligned.out.bam  -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule star_pe:
    input:
        trim_fq1="{output_dir}/{{sample}}/trim/{trim_soft}/{{sample}}_1_val_1.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
        trim_fq2="{output_dir}/{{sample}}/trim/{trim_soft}/{{sample}}_2_val_2.fq.gz".format(output_dir=config['output_dir'],trim_soft=config['trim_soft']),
    output:
        bam="{output_dir}/{{sample}}/bam/star/{{sample}}_pair.bam".format(output_dir=config['output_dir']),
    log:
        "{output_dir}/{{sample}}/logs/star/{{sample}}.log".format(output_dir=config['output_dir']),
    threads: 10
    conda:
        "flask2024",
    params:
        ref=config["ref_path"] +"/"+config["ref_species"] +"/"+ config["ref_release"] +"/star/"+config["ref_release"]+"_"+config["ref_len"],
        bam="{output_dir}/{{sample}}/bam/{{sample}}.star.".format(output_dir=config['output_dir']),
        other=config['align_params'],
    shell:
        """
        ({STAR} --genomeDir {params.ref} --readFilesIn {input.trim_fq1} {input.trim_fq2}  --outFileNamePrefix {params.bam} --runThreadN {threads} \
             --outSAMtype BAM Unsorted --outSAMmultNmax 1  --readFilesCommand zcat {params.other} ) > {log} 2>&1
        samtools sort {params.bam}Aligned.out.bam -@ {threads} -o {output.bam}
        samtools index {output.bam}


        """
