import argparse
import collections
import gzip
import multiprocessing
import subprocess
from collections import Counter
import pysam
from Bio import SeqIO
from multiprocessing import Pool


class BAMSplitter:
    @classmethod
    def get_total_lines(cls, bam_file):
        result = subprocess.run(
            f"samtools view {bam_file} | wc -l",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"Failed to get total lines: {result.stderr}")
        return int(result.stdout.strip())

    @classmethod
    def split_bam(cls, params):
        bam_file, start, end, chunk_index = params
        output_bam = f"{bam_file[:-4]}_chunk_{chunk_index}.bam"
        command = (
            f"samtools view {bam_file} | "
            f"sed -n '{start},{end}p' | "
            f"cat {bam_file}.header.sam - | "
            f"samtools view -Sb - > {output_bam}"
        )
        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"samtools index {output_bam}", shell=True, check=True)
        print(f"Generated and indexed {output_bam}")

    @classmethod
    def split_bam_file(cls, bam_file, num_chunks):
        total_lines = cls.get_total_lines(bam_file)
        lines_per_chunk = (total_lines + num_chunks - 1) // num_chunks
        subprocess.run(f"samtools view -H {bam_file} > {bam_file}.header.sam", shell=True, check=True)
        split_params = []
        for i in range(num_chunks):
            start = i * lines_per_chunk + 1
            end = min((i + 1) * lines_per_chunk, total_lines)
            split_params.append((bam_file, start, end, i + 1))
        with Pool(processes=num_chunks) as pool:
            pool.map(cls.split_bam, split_params)
        print("All chunks have been processed and indexed.")


class DNASeqUtility:
    @staticmethod
    def readgenome(name):
        dd = {}
        with open(name, 'r') as input_fasta:
            for record in SeqIO.parse(input_fasta, 'fasta'):
                dd[record.id] = str(record.seq).upper()
        return dd

    @staticmethod
    def reverse_complement(dna):
        complement = str.maketrans('ATCG', 'TAGC')
        reverse = dna[::-1].translate(complement)
        return reverse

    @staticmethod
    def get_seq_dir(fastq_file):
        seq_dict = {}
        opener = gzip.open if fastq_file.endswith('.gz') else open
        with opener(fastq_file, 'rt') as f:
            while True:
                line = f.readline().strip()
                if not line:
                    break
                if line.startswith('@'):
                    seq_id = line.split()[0][1:]
                    seq = f.readline().strip()
                    f.readline()
                    f.readline()
                    seq_dict[seq_id] = seq
        return seq_dict

    @staticmethod
    def get_complementary_base(base):
        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c"
        }
        return complement.get(base, base)

    @staticmethod
    def safe_get_base(ref_dict, reference_name, ref_pos, is_reverse):
        if ref_pos is None or ref_pos < 0 or ref_pos >= len(ref_dict[reference_name]):
            return None
        elif is_reverse:
            return DNASeqUtility.get_complementary_base(ref_dict[reference_name][ref_pos])
        else:
            return ref_dict[reference_name][ref_pos]

    @staticmethod
    def process_aligned_pairs(read, ref_dict):
        aligned_pairs = read.get_aligned_pairs(matches_only=False)
        is_reverse = read.is_reverse
        if read.is_read2:
            if is_reverse:
                align_list = [ref_pos for query_pos, ref_pos in aligned_pairs if query_pos is not None]
                last_pos = next((ref_pos for ref_pos in reversed(align_list) if ref_pos is not None), None)
                if last_pos is None:
                    return []
                align_list.extend([last_pos + 1, last_pos + 2])
                result = [DNASeqUtility.safe_get_base(ref_dict, read.reference_name, ref_pos, False) for ref_pos in
                          align_list]
            else:
                align_list = [ref_pos for query_pos, ref_pos in aligned_pairs if query_pos is not None][::-1]
                last_pos = next((ref_pos for ref_pos in reversed(align_list) if ref_pos is not None),None)
                if last_pos is None:
                    return []
                align_list.extend([last_pos - 1, last_pos - 2])
                result = [DNASeqUtility.safe_get_base(ref_dict, read.reference_name, ref_pos, True) for ref_pos in
                          align_list]
        else:
            if is_reverse:
                align_list = [ref_pos for query_pos, ref_pos in aligned_pairs if query_pos is not None][::-1]
                last_pos = next((ref_pos for ref_pos in reversed(align_list) if ref_pos is not None),None)
                if last_pos is None:
                    return []
                align_list.extend([last_pos - 1, last_pos - 2])
                result = [DNASeqUtility.safe_get_base(ref_dict, read.reference_name, ref_pos, True) for ref_pos in
                          align_list]
            else:
                align_list = [ref_pos for query_pos, ref_pos in aligned_pairs if query_pos is not None]
                last_pos = next((ref_pos for ref_pos in reversed(align_list) if ref_pos is not None),None)
                if last_pos is None:
                    return []
                align_list.extend([last_pos + 1, last_pos + 2])
                result = [DNASeqUtility.safe_get_base(ref_dict, read.reference_name, ref_pos, False) for ref_pos in
                          align_list]
        return result


# 基础规则类
class BaseRule:
    def __init__(self, name):
        self.name = name  # 规则名称

    def apply(self, read1, read2, align_next):
        """子类必须实现的规则应用方法"""
        raise NotImplementedError("必须在子类中实现 apply 方法")

    def get_XM_char(self, ref_list, xm_values):
        if xm_values == ".":
            return "."
        if ref_list[0] is None:
            return xm_values[3]
        elif ref_list[0].upper() == "G":
            return xm_values[0]
        elif ref_list[1] is None:
            return xm_values[3]
        elif ref_list[0].upper() == "N" or ref_list[1].upper() == "N":
            return xm_values[3]
        elif ref_list[1].upper() == "G":
            return xm_values[1]
        return xm_values[2]
class MatrixRule1(BaseRule):
    def __init__(self):
        super().__init__("MatrixRule1")
        self.base_to_index = {"A": 0, "G": 1, "C": 2, "T": 3, "N": 4}
        self.rule_matrix = [
            # A    G    C    T    N
            [".", ".", ".", ".", "."],  # A
            [".", ".", ".", ".", "."],  # G
            [".", ".", ["z", "x", "h", "u"], ["y", "v", "w", "n"], ["y", "v", "w", "n"]],  # C
            [".", ".", ["Z", "X", "H", "U"], ["Y", "V", "W", "N"], ["y", "v", "w", "n"]],  # T
            [".", ".", ["y", "v", "w", "n"], ["y", "v", "w", "n"], ["y", "v", "w", "n"]]  # N
        ]

    def apply(self, read1, read2, align_next):
        idx1 = self.base_to_index.get(read1, 4)
        idx2 = self.base_to_index.get(read2, 4)
        choices = self.rule_matrix[idx1][idx2]
        return self.get_XM_char(align_next, choices)


# 第二种规则：矩阵规则2
class MatrixRule2(BaseRule):
    # 这里事read为T的作为修饰，也就是基于转化的
    def __init__(self):
        super().__init__("MatrixRule2")
        self.base_to_index = {"A": 0, "G": 1, "C": 2, "T": 3, "N": 4}
        self.rule_matrix = [
            # A    G    C    T    N
            [".", ".", ".", ".", "."],  # A
            [".", ".", ".", ".", "."],  # G
            [".", ".", ["z", "x", "h", "u"], ".", "."],  # C
            [".", ".", ".", ["Z", "X", "H", "U"], "."],  # T
            [".", ".", ".", ".", "."]  # N
        ]

    def apply(self, read1, read2, align_next):
        idx1 = self.base_to_index.get(read1, 4)
        idx2 = self.base_to_index.get(read2, 4)
        choices = self.rule_matrix[idx1][idx2]
        return self.get_XM_char(align_next, choices)
class RuleFactory:
    @staticmethod
    def get_rule(rule_key):
        if rule_key in {"CGS"}:
            return MatrixRule1()
        elif rule_key in {"CT"}:
            return MatrixRule2()
        else:
            raise ValueError(f"未定义的规则: {rule_key}")







def call_bam_xm(bam_file):
    bed_file = bam_file[:-4] + ".bed"
    update_map = {
        "Z": (0, 1, 0, 0), "z": (0, 0, 0, 1), "X": (0, 1, 0, 0), "x": (0, 0, 0, 1),
        "H": (0, 1, 0, 0), "h": (0, 0, 0, 1), "U": (0, 1, 0, 0), "u": (0, 0, 0, 1),
        "Y": (1, 0, 0, 0), "y": (0, 0, 1, 0), "V": (1, 0, 0, 0), "v": (0, 0, 1, 0),
        "W": (1, 0, 0, 0), "w": (0, 0, 1, 0), "N": (1, 0, 0, 0), "n": (0, 0, 1, 0),
    }
    bed_dict = collections.defaultdict(lambda: (0, 0, 0, 0))
    bam = pysam.AlignmentFile(bam_file, "r")
    out_file = open(bed_file, "w")
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        chrom = read.reference_name
        xm_value = read.get_tag("XM")
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

        for query_pos, ref_pos in aligned_pairs:
            if query_pos is None or ref_pos is None:
                continue
            if xm_value[query_pos] in update_map:
                context = (
                    "CpG" if xm_value[query_pos] in "ZzYy" else
                    "CHG" if xm_value[query_pos] in "XxVv" else
                    "CHH" if xm_value[query_pos] in "HhWw" else
                    "CN"
                )
                key = (chrom, ref_pos, context)
                bed_dict[key] = tuple(map(sum, zip(bed_dict[key], update_map[xm_value[query_pos]])))

    for key, counts in bed_dict.items():
        chrom, ref_pos, context = key
        if ref_pos is None:  # 额外的保护措施
            continue
        hmc_counts = counts[0]
        mc_counts = counts[1]
        unhmc_counts = counts[2]
        C_counts = counts[3]
        total = counts[0] + counts[1] + counts[3]
        if total == 0:
            out_file.write(
                f"{chrom}\t{ref_pos + 1}\t{ref_pos + 1}\t{0:.4f}\t{0:.4f}\t{hmc_counts}\t{mc_counts}\t{unhmc_counts}\t{C_counts}\t{context}\n"
            )
        else:
            hmc_ratio = counts[0] / total
            mc_ratio = counts[1] / total
            out_file.write(
                f"{chrom}\t{ref_pos + 1}\t{ref_pos + 1}\t{hmc_ratio:.4f}\t{mc_ratio:.4f}\t{hmc_counts}\t{mc_counts}\t{unhmc_counts}\t{C_counts}\t{context}\n"
            )
    bam.close()
    out_file.close()


def extract_XM_modify(bam_file, out_file, ref_dict, read1_dict, read2_dict, rule_key):
    output_bam_file = out_file
    with pysam.AlignmentFile(bam_file, "r") as bam, pysam.AlignmentFile(output_bam_file, "wb",
                                                                        header=bam.header) as out_bam:
        rule = RuleFactory.get_rule(rule_key)
        print(rule.rule_matrix)
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            align_list = DNASeqUtility.process_aligned_pairs(read, ref_dict)
            # 这里有可能返回空

            if not align_list:  # 当 align_list 为空时
                read.set_tag("XM", "")  # 设置 XM 为空字符串

            elif not read1_dict and not read2_dict:

                read_seq = read.query_sequence
                if read.is_read2:
                    read_seq = DNASeqUtility.reverse_complement(read_seq)
                if read.is_reverse:
                    read_seq = DNASeqUtility.reverse_complement(read_seq)

                XM = [
                    rule.apply(read_seq[i], read_seq[i], [align_list[i + 1], align_list[i + 2]])
                    if align_list[i] == "C"
                    else "."
                    for i in range(len(align_list) - 2)
                ]
                xm_value = XM[::-1] if read.is_reverse else XM
                read.set_tag("XM", "".join(xm_value[::-1] if read.is_read2 else xm_value))
            else:
                read1 = read1_dict[read.query_name]
                read2 = read2_dict[read.query_name]
                XM = [
                    rule.apply(read1[i], read2[i], [align_list[i + 1], align_list[i + 2]])
                    if align_list[i] == "C"
                    else "."
                    for i in range(len(align_list) - 2)
                ]
                read.set_tag("XM", "".join(XM[::-1] if read.is_reverse else XM))
            out_bam.write(read)
    subprocess.run(f"samtools index {output_bam_file}", shell=True, check=True)
    call_bam_xm(output_bam_file)


def process_bam_wrapper(params):
    bam_file, ref_dict, read1_dict, read2_dict, rule_key = params
    extract_XM_modify(bam_file, ref_dict, read1_dict, read2_dict, rule_key)
def batch_process_bam(bam_files, ref_dict, read1_dict, read2_dict, rule_key, num_processes):
    params = [(bam_file, ref_dict, read1_dict, read2_dict, rule_key) for bam_file in bam_files]
    with Pool(processes=num_processes) as pool:
        pool.map(process_bam_wrapper, params)








def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--parallel",
        help="Number of parallel threads to use for processing. Default: 20",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--ref",
        help="ref fasta",
        default="/data/wangzc/cgs/mc/new_mm9_lam_mc.fa",
        type=str,
    )
    parser.add_argument(
        "--bam",
        help="bam file",
        type=str,
    )
    parser.add_argument(
        "--output",
        help="output file",
        type=str,
    )
    parser.add_argument(
        "--cutfq1",
        help="output file",
        type=str,
        default="",
    )
    parser.add_argument(
        "--cutfq2",
        help="output file",
        type=str,
        default="",
    )
    parser.add_argument(
        "--rule",
        help="Rule key to specify which processing rule to use (e.g., CT, CGS)",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    if args.cutfq1 == "":
        read1_dict = {}
    else:
        read1_dict = DNASeqUtility.get_seq_dir(args.cutfq1)
    if args.cutfq2 == "":
        read2_dict = {}
    else:
        read2_dict = DNASeqUtility.get_seq_dir(args.cutfq2)
    ref_dict = DNASeqUtility.readgenome(args.ref)



    if args.parallel == 1:
        extract_XM_modify(args.bam, args.output,ref_dict, read1_dict, read2_dict, args.rule)

    else:
        BAMSplitter.split_bam_file(args.bam, args.parallel)
        bam_files = [f"{args.bam[:-4]}_chunk_{i}.bam" for i in range(1, args.parallel + 1)]
        batch_process_bam(bam_files, ref_dict, read1_dict, read2_dict, args.rule, args.parallel)



if __name__ == '__main__':
    main()











