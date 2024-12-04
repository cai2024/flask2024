import argparse
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
    if cover * ref_size == 0:
        average_depth = 0
    else:
        average_depth = round(map_size / (cover * ref_size * 0.01), 4)
    return cover, average_depth, round(map_size / (ref_size ), 4)


def sample_sort_and_index_bam(input_bam, output_file, num, qualimap="qualimap"):
    result_cover = []
    result_depth = []
    result_x = []
    result_file = f"{output_file}.txt"
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))  # 创建输出目录，如果它不存在
    for i in range(1, num + 1):
        output_bam = f"{output_file}_sampled_sorted_{i}_of_{num}.bam"
        output_cover = f"{output_file}_sampled_sorted_{i}_of_{num}"
        cover_file = os.path.join(output_cover, "genome_results.txt")
        if i == num:
            command_sort = f"samtools view -b {input_bam} | samtools sort -@ 10 -o {output_bam} && samtools index -@ 10 {output_bam}"
        else:
            command_sort = f"samtools view -s 38.{str(round(i / num, 4))[2:]} -b {input_bam} | samtools sort -@ 10 -o {output_bam} && samtools index -@ 10 {output_bam}"
        subprocess.run(command_sort, shell=True, check=True)
        command_qualimap = [qualimap, "bamqc", "-bam", output_bam, "-outdir", output_cover, "-outformat", "PDF:HTML", "-nt", "10", "--java-mem-size=10G"]
        subprocess.run(command_qualimap, check=True)
        c, d, e = cover_depth(cover_file)
        result_cover.append(str(c))
        result_depth.append(str(d))
        result_x.append(str(e))
    with open(result_file, "w") as file:
        print(' '.join(result_cover), file=file)
        print(' '.join(result_depth), file=file)
        print(' '.join(result_x), file=file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        help="bam file",
        type=str,
    )
    parser.add_argument(
        "--output",
        help="output file",
        default="N",
        type=str,
    )
    parser.add_argument(
        "--num",
        help="num",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--qualimap",
        help="qualimap",
        default="qualimap",
        type=str,
    )
    args = parser.parse_args()
    sample_sort_and_index_bam(args.bam, args.output, args.num, args.qualimap)


