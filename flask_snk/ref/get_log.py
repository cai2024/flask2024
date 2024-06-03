import os
import subprocess
import json
import argparse


def spikein_log(spikein, bam, hmc="" , mc="", bedtools="bedtools"):
    bed = bam[:-3] + "bed"
    command = f"{bedtools} genomecov -ibam {bam} -bga > {bed}"
    process = subprocess.run(command, shell=True, text=True, check=True)

    spikein = spikein.split(',')
    ref={"lambda":48502,
         "mc":897,
         "G5hmc":335,
         "puc19":2686,
    }
    cover_base = {key: 0 for key in spikein}
    depth_base = {key: 0 for key in spikein}


    with open(bed) as file:
        for line in file:
            line = line.split()
            if int(float(line[3])) != 0 and line[0] in spikein:
                cover_base[line[0]] += (int(line[2]) - int(line[1]))
                depth_base[line[0]] += ((int(line[2]) - int(line[1])) * int(float(line[3])))
    result=[]
    if hmc=="":
        for key in spikein:
            if key in ref:
                cover=str(round(cover_base[key] / ref[key], 4))
                depth="0" if cover_base[key] == 0 else str(round(depth_base[key] / cover_base[key], 4))
                result.append(f'{cover}|{depth}')
        print('\t'.join(result))

    else:
        for key in spikein:
            if key in ref:
                cover=str(round(cover_base[key] / ref[key], 4))
                depth="0" if cover_base[key] == 0 else str(round(depth_base[key] / cover_base[key], 4))
                command1 = f"""awk '$7 == "CpG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {hmc}"""
                command2 = f"""awk '$7 == "CHG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {hmc}"""
                command3 = f"""awk '$7 == "CHH" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {hmc}"""
                temp1 = subprocess.run(command1, shell=True, text=True, capture_output=True)
                temp2 = subprocess.run(command2, shell=True, text=True, capture_output=True)
                temp3 = subprocess.run(command3, shell=True, text=True, capture_output=True)
                cpg = 0 if temp1.stdout.strip()=="" else round(float(temp1.stdout.strip()),4)
                chg = 0 if temp2.stdout.strip() == "" else round(float(temp2.stdout.strip()), 4)
                chh = 0 if temp3.stdout.strip() == "" else round(float(temp3.stdout.strip()), 4)
                command1 = f"""awk '$7 == "CpG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {mc}"""
                command2 = f"""awk '$7 == "CHG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {mc}"""
                command3 = f"""awk '$7 == "CHH" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {mc}"""
                temp1 = subprocess.run(command1, shell=True, text=True, capture_output=True)
                temp2 = subprocess.run(command2, shell=True, text=True, capture_output=True)
                temp3 = subprocess.run(command3, shell=True, text=True, capture_output=True)
                mc_cpg = 0 if temp1.stdout.strip()=="" else round(float(temp1.stdout.strip()),4)
                mc_chg = 0 if temp2.stdout.strip() == "" else round(float(temp2.stdout.strip()), 4)
                mc_chh = 0 if temp3.stdout.strip() == "" else round(float(temp3.stdout.strip()), 4)
                result.append(f'{cover}|{depth}\t{cpg}|{chg}|{chh}\t{mc_cpg}|{mc_chg}|{mc_chh}')
        print('\t'.join(result))


def spikein_wgbs(spikein, bam , sort_bam ,mc="", bedtools="bedtools"):
    bed = bam[:-3] + "bed"
    command = f"{bedtools} genomecov -ibam {sort_bam} -bga > {bed}"
    process = subprocess.run(command, shell=True, text=True, check=True)

    spikein = spikein.split(',')
    ref={"lambda":48502,
         "mc":897,
         "G5hmc":335,
         "puc19":2686,
    }
    cover_base = {key: 0 for key in spikein}
    depth_base = {key: 0 for key in spikein}


    with open(bed) as file:
        for line in file:
            line = line.split()
            if int(float(line[3])) != 0 and line[0] in spikein:
                cover_base[line[0]] += (int(line[2]) - int(line[1]))
                depth_base[line[0]] += ((int(line[2]) - int(line[1])) * int(float(line[3])))
    result=[]
    if mc=="":
        for key in spikein:
            if key in ref:
                cover=str(round(cover_base[key] / ref[key], 4))
                depth="0" if cover_base[key] == 0 else str(round(depth_base[key] / cover_base[key], 4))
                result.append(f'{cover}|{depth}')
        print('\t'.join(result))

    else:
        for key in spikein:
            if key in ref:
                cover=str(round(cover_base[key] / ref[key], 4))
                depth="0" if cover_base[key] == 0 else str(round(depth_base[key] / cover_base[key], 4))
                command1 = f"""awk '$6 == "CG" && $1 == "{key}" && ($4 + $5) != 0 {{ sum += $4 / ($4 + $5); count++ }} END {{ if (count > 0) print sum / count; else print 0;}}' {mc}"""
                command2 = f"""awk '$6 == "CHG" && $1 == "{key}" && ($4 + $5) != 0 {{ sum += $4 / ($4 + $5); count++ }} END {{ if (count > 0) print sum / count; else print 0;}}' {mc}"""
                command3 = f"""awk '$6 == "CHH" && $1 == "{key}" && ($4 + $5) != 0 {{ sum += $4 / ($4 + $5); count++ }} END {{ if (count > 0) print sum / count; else print 0;}}' {mc}"""
                temp1 = subprocess.run(command1, shell=True, text=True, capture_output=True)
                temp2 = subprocess.run(command2, shell=True, text=True, capture_output=True)
                temp3 = subprocess.run(command3, shell=True, text=True, capture_output=True)
                mc_cpg = 0 if temp1.stdout.strip()=="" else round(float(temp1.stdout.strip()),4)
                mc_chg = 0 if temp2.stdout.strip() == "" else round(float(temp2.stdout.strip()), 4)
                mc_chh = 0 if temp3.stdout.strip() == "" else round(float(temp3.stdout.strip()), 4)
                result.append(f'{cover}|{depth}\t{mc_cpg}|{mc_chg}|{mc_chh}')
        print('\t'.join(result))


def size_clean(raw_fq1, trim_fq1, mode):

    if mode == "single":
        raw_fq1 = raw_fq1 + ".fq.gz"
        trim_fq1 = trim_fq1 + "_trimmed.fq.gz"
    else:
        raw_fq1 = raw_fq1 + "_1.fq.gz"
        trim_fq1 = trim_fq1 + "_1_val_1.fq.gz"


    command = f"seqtk comp {raw_fq1} | awk '{{sum += $2}} END {{print sum}}'"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if mode == "single":
        sequence_size = round(int(float(result.stdout.strip())) / 1000000000, 4)
    else:
        sequence_size = round(2 * int(float(result.stdout.strip())) / 1000000000, 4)

    command = f"seqtk seq {raw_fq1} | wc -l"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    raw_read_num = int(result.stdout.strip()) // 4
    command = f"seqtk seq {trim_fq1} | wc -l"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    clean_read_num = int(result.stdout.strip()) // 4
    clean_read_ratio = round(clean_read_num / raw_read_num, 4)
    print(f"{sequence_size}\t{clean_read_ratio}")


def get_json(json_file, json_key):
    if os.path.exists(json_file):
        key_list = json_key.split(',')
        with open(json_file, 'r') as file:
            data = json.load(file)
            if len(key_list) == 2:
                my_data = data[key_list[0]][key_list[1]]
            elif len(key_list) == 3:
                my_data = data[key_list[0]][key_list[1]][key_list[2]]
        print(my_data)
    else:
        print("0")



def dup_map(fastq, mode, bam, debam, dup_radio1):
    command = f"samtools view -c -F 4 {bam}"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    map_read_num=int(result.stdout.strip())

    command = f"samtools view -c -F 4 {debam}"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    de_read_num=int(result.stdout.strip())

    command = f"zcat {fastq} | wc -l"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    all_read_num=int(result.stdout.strip())/2 if mode == "pair" else int(result.stdout.strip())/4
    
    if map_read_num==0:
        duplication_rate = round(float(dup_radio1), 4)
    else:
        d1 = round(1 - round(de_read_num/map_read_num, 4),4)
        duplication_rate=  str(round(float(dup_radio1),4))+ "|" + str(d1)
    if all_read_num==0:
        map_read_ratio=0
    else:
        map_read_ratio=round(map_read_num/all_read_num,4)
    print(f"{duplication_rate}\t{map_read_ratio}")

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

    print(f"{cover}\t{average_depth}")







def mit_ratio(bam):
    command = f"samtools view -c {bam} chrM"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    mit_read_num=int(result.stdout.strip())

    command = f"samtools view -c -F 4 {bam}"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    map_read_num=int(result.stdout.strip())

    if map_read_num==0:
        print("0")
    else:
        print(round(mit_read_num/map_read_num, 4))



def gc_q20_q30(fq_path, mode, bam, picard, java = "java"):
    if mode=="single":
        fq1 = fq_path + "_trimmed.fq.gz"
        command = f"seqtk fqchk {fq1} | awk 'NR == 3 {{print $4 + $5}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        gc = round(float(result.stdout.strip()),4)
    elif mode=="pair":
        fq1 = fq_path + "_1_val_1.fq.gz"
        fq2 = fq_path + "_2_val_2.fq.gz"
        command = f"seqtk fqchk {fq1} | awk 'NR == 3 {{print $4 + $5}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        gc1 = float(result.stdout.strip())
        command = f"seqtk fqchk {fq2} | awk 'NR == 3 {{print $4 + $5}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        gc2 = float(result.stdout.strip())
        gc = round((gc1+gc2)/2,4)
    elif mode=="CGS":
        fq1 = fq_path
        command = f"seqtk fqchk {fq1} | awk 'NR == 3 {{print $4 + $5}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        gc = round(float(result.stdout.strip()),4)


    command = f"{java} -jar {picard} CollectQualityYieldMetrics -I {bam} -O /dev/stdout"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    my_list = str(result).split('\\n')
    indexes = [index for index, s in enumerate(my_list) if s.startswith("TOTAL_READS")]
    qc_list = my_list[indexes[0]+1].split("\\t")

    q20=round(int(qc_list[5])/int(qc_list[3]),4)
    q30=round(int(qc_list[7])/int(qc_list[3]),4)

    print(f"{gc}\t{q20}\t{q30}")







def main():
    parser = argparse.ArgumentParser(description='Run different bioinformatics tools.')
    subparsers = parser.add_subparsers(dest='command')

    # spikein 子命令
    parser_spikein = subparsers.add_parser('spikein', help='Run spikein')
    parser_spikein.add_argument('--spikein', type=str, help='The input file for FastQC.')
    parser_spikein.add_argument('--bam', type=str, help='bam.')
    parser_spikein.add_argument('--hmc', type=str, default='', help='Optional: Path to the HMC file.')
    parser_spikein.add_argument('--mc', type=str, default='', help='Optional: Path to the MC file.')
    parser_spikein.add_argument('--bedtools', type=str, default='bedtools', help='Optional: Path to the bedtools executable or name if in PATH.')
    # spikein_wgbs 子命令
    parser_spikein_wgbs = subparsers.add_parser('spikein_wgbs', help='Run spikein')
    parser_spikein_wgbs.add_argument('--spikein', type=str, help='The input file for FastQC.')
    parser_spikein_wgbs.add_argument('--bam', type=str, help='bam.')
    parser_spikein_wgbs.add_argument('--sort_bam', type=str, help='bam.')
    parser_spikein_wgbs.add_argument('--mc', type=str, default='', help='Optional: Path to the HMC file.')
    parser_spikein_wgbs.add_argument('--bedtools', type=str, default='bedtools', help='Optional: Path to the bedtools executable or name if in PATH.')

    # size_clean 子命令
    parser_size_clean = subparsers.add_parser('size_clean', help='Run size_clean')
    parser_size_clean.add_argument('--raw_fq1', type=str, help='The input ')
    parser_size_clean.add_argument('--trim_fq1', type=str, help='The input ')
    parser_size_clean.add_argument('--mode', type=str, help='The input ')

    # get_json 子命令
    parser_get_json = subparsers.add_parser('get_json', help='Run size_clean')
    parser_get_json.add_argument('--json_file', type=str, help='The input ')
    parser_get_json.add_argument('--json_key', type=str, help='The input ')

    # dup_map 子命令
    parser_dup_map = subparsers.add_parser('dup_map', help='Run size_clean')
    parser_dup_map.add_argument('--fastq', type=str, help='The input ')
    parser_dup_map.add_argument('--mode', type=str, help='The input ')
    parser_dup_map.add_argument('--bam', type=str, help='The input ')
    parser_dup_map.add_argument('--debam', type=str, help='The input ')
    parser_dup_map.add_argument('--dup_radio1', type=str, help='The input ')

    # cover_depth 子命令
    parser_cover_depth = subparsers.add_parser('cover_depth', help='Run size_clean')
    parser_cover_depth.add_argument('--input_file', type=str, help='The input ')

    # mit_ratio 子命令
    parser_mit_ratio = subparsers.add_parser('mit_ratio', help='Run size_clean')
    parser_mit_ratio.add_argument('--bam', type=str, help='The input ')

    # gc_q20_q30 子命令
    parser_gc_q20_q30 = subparsers.add_parser('gc_q20_q30', help='Run size_clean')
    parser_gc_q20_q30.add_argument('--fq_path', type=str, help='The input ')
    parser_gc_q20_q30.add_argument('--mode', type=str, help='The input ')
    parser_gc_q20_q30.add_argument('--bam', type=str, help='The input ')
    parser_gc_q20_q30.add_argument('--picard', type=str, help='The input ')
    parser_gc_q20_q30.add_argument('--java', type=str, help='The input ')

    args = parser.parse_args()



    if args.command == 'spikein':
        spikein_log(args.spikein, args.bam, args.hmc, args.mc, args.bedtools)
    elif args.command == 'spikein_wgbs':
        spikein_wgbs(args.spikein, args.bam, args.sort_bam ,args.mc, args.bedtools)
    elif args.command == 'size_clean':
        size_clean(args.raw_fq1, args.trim_fq1, args.mode)
    elif args.command == 'get_json':
        get_json(args.json_file, args.json_key)
    elif args.command == 'dup_map':
        dup_map(args.fastq, args.mode, args.bam, args.debam, args.dup_radio1)
    elif args.command == 'cover_depth':
        cover_depth(args.input_file)
    elif args.command == 'mit_ratio':
        mit_ratio(args.bam)
    elif args.command == 'gc_q20_q30':
        gc_q20_q30(args.fq_path,args.mode,args.bam,args.picard,args.java)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

