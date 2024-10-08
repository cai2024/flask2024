import subprocess
import os

def cover_depth(input_file):
    command = f"sed -n 's/.*There is a \\([0-9.]*\\)% of reference with a coverageData >= 1X.*/\\1/p' {input_file}"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    cover = float(result.stdout.strip())


    command = f"sed -n 's/.*number of bases = \\([0-9,]*\\) bp.*/\\1/p' {input_file} | tr -d ','"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    ref_size = int(result.stdout.strip())

    command = f"sed -n 's/.*number of mapped bases = \\([0-9,]*\\) bp.*/\\1/p' {input_file} | tr -d ','"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    map_size = int(result.stdout.strip())

    if cover * ref_size ==0:
        average_depth = 0
    else:
        average_depth = round(map_size/(cover * ref_size * 0.01), 4)
    
    return cover,average_depth


def sample_sort_and_index_bam(input_bam, n, output_dir, output_prefix):
    """
    对 BAM 文件进行抽样、排序并建立索引，使用多线程来加速处理，并指定输出文件夹位置。
    
    :param input_bam: 输入 BAM 文件路径
    :param n: 抽样的分母整数
    :param output_dir: 输出文件夹路径
    :param output_prefix: 输出文件的前缀名
    """
    result_cover=[]
    result_depth=[]
    result_file = os.path.join(output_dir, f"cover_depth.txt")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # 创建输出目录，如果它不存在
    
    for i in range(1, n+1):
        output_bam = os.path.join(output_dir, f"{output_prefix}_sampled_sorted_{i}_of_{n}.bam")
        output_cover = os.path.join(output_dir, f"{output_prefix}_sampled_sorted_{i}_of_{n}")
        cover_file = os.path.join(output_cover,"genome_results.txt")
        if i == n:
            # 最后一个分数，复制原始文件并排序
            command_sort = f"samtools view -b {input_bam} | samtools sort -@ 10 -o {output_bam}"
        else:
            # 计算抽样概率并抽样排序
            fraction = i / n
            command_sort = f"samtools view -s 42.{int(fraction*100)} -b {input_bam} | samtools sort -@ 10 -o {output_bam}"
        
        # 执行排序命令
        print(f"Running: {command_sort}")
        subprocess.run(command_sort, shell=True, check=True)
        print(f"Sorted: {output_bam}")

        # 执行索引命令
        command_index = f"samtools index -@ 10 {output_bam}"
        print(f"Running: {command_index}")
        subprocess.run(command_index, shell=True, check=True)
        print(f"Indexed: {output_bam}.bai")

        # 运行 qualimap bamqc 命令
        command_qualimap = [
            "/data/wangzc/soft/qualimap_v2.3/qualimap", "bamqc", "-bam", output_bam, "-outdir", output_cover,
            "-outformat", "PDF:HTML", "-nt", "10", "--java-mem-size=10G"
        ]
        print(f"Running: {' '.join(command_qualimap)}")
        subprocess.run(command_qualimap, check=True)
        
        
        c,d = cover_depth(cover_file)
        result_cover.append(c)
        result_depth.append(d)

    result_cover = [str(item) for item in result_cover]  
    result_depth = [str(item) for item in result_depth]   
    with open(result_file,"w")as file:
        print( ' '.join(result_cover),file=file)
        print( ' '.join(result_depth),file=file)




# 示例调用
input_bam = "CGS-13.5_dedup.bam"
n = 5  # 用户提供的整数
output_dir = "/data/all_data/0628/CGS-13.5/data/dedup/test"  # 指定输出文件夹
output_prefix = "mydata"
sample_sort_and_index_bam(input_bam, n, output_dir, output_prefix)

