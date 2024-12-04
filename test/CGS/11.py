import Levenshtein
import gzip
import math
from functools import partial
import multiprocessing
from Bio.Align.substitution_matrices import Array
from Bio.Align import PairwiseAligner
import numpy as np
import argparse
from abc import ABC, abstractmethod
from Bio import SeqIO
from multiprocessing import Pool


class DNASeqUtility:
    @staticmethod
    def read_genome(name):
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


class FileReader:
    def __init__(self, file1, file2):
        # 判断文件是否是 .gz 压缩文件
        self.file1 = self._open_file(file1)
        self.file2 = self._open_file(file2)

    def _open_file(self, file_path):
        if file_path.endswith(".gz"):
            return gzip.open(file_path, "rt")
        else:
            return open(file_path, "r")

    def __iter__(self):
        return self

    def __next__(self):
        # 每次读取 4 行（FASTQ 格式）
        lines1 = [self.file1.readline().strip() for _ in range(4)]
        lines2 = [self.file2.readline().strip() for _ in range(4)]
        # 检查文件结束情况
        if not lines1[0] or not lines2[0]:
            raise StopIteration
        # 返回 PairedRead 对象
        return PairedRead(
            name=lines1[0].split()[0],
            read1=lines1[1],
            read2=lines2[1],
            qual1=lines1[3],
            qual2=lines2[3],
        )

    def close(self):
        # 关闭文件句柄
        self.file1.close()
        self.file2.close()




class TrimStrategy(ABC):
    @staticmethod
    def check(str1, str2):
        """
        str1‘end == str2’start
        xxxxxxxxxxxxxAAA
        AAAxxxxxxxxxxxxx
        :return 3
        """
        length = min(len(str1), len(str2))
        k = max(range(length + 1),
                key=lambda i: i if Levenshtein.hamming(str1[-i:], str2[:i]) < i * 0.1 else False)
        return k

    @abstractmethod
    def trim(self, read1, read2):
        """定义trim_overlap的接口"""
        return read1, read2

class DefaultTrimStrategy(TrimStrategy):
    def trim(self, seq1, seq2):
        reverse_seq2 = DNASeqUtility.reverse_complement(seq2)
        overlap_len1 = self.check(seq1, reverse_seq2)
        if overlap_len1 > 5:
            trim_length = int(overlap_len1 / 2)
            return seq1[:-trim_length], seq2[:-trim_length]

        overlap_len2 = self.check(reverse_seq2, seq1)
        if overlap_len2 > min(len(seq1), len(seq2)) * 0.8 and overlap_len2 > 5:
            trim_length = int(overlap_len2 / 2)
            return seq1[:-trim_length], seq2[:-trim_length]

        return seq1, seq2

class AdvancedTrimStrategy(TrimStrategy):
    def trim(self, seq1, seq2):
        # 这里业务逻辑未确定，临时占位
        overlap_len1 = self.check(seq1, DNASeqUtility.reverse_complement(seq2))
        if overlap_len1 > 5:
            trim_length = int(overlap_len1 / 2)
            return seq1[:-trim_length], seq2[:-trim_length]

class TrimStrategyFactory:
    @staticmethod
    def get_strategy(strategy_type):
        if strategy_type == "default":
            return DefaultTrimStrategy()
        elif strategy_type == "advanced":
            return AdvancedTrimStrategy()
        else:
            raise ValueError("Unknown strategy type")






class AlignmentConfig:
    """提供不同模式的替换矩阵和打分规则"""

    @staticmethod
    def get_rule_matrix_ct():
        return np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0]
        ])

    @staticmethod
    def get_rule_matrix_ctga():
        return np.array([
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0]
        ])

    @staticmethod
    def get_score_rules_ct():
        return {
            ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
            ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
            ('G', 'A'): -20, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
            ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
            ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
            ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
        }

    @staticmethod
    def get_score_rules_ctga():
        return {
            ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
            ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
            ('G', 'A'): 1, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
            ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
            ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
            ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
        }

class BaseAlignment(ABC):
    def __init__(self, rule_matrix, score_rules, mode='global', open_gap_score=-3, extend_gap_score=-2):
        self.alphabet = "ACGTN"
        self.aligner = PairwiseAligner()
        self.aligner.mode = mode
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
        self.rule_matrix = rule_matrix
        self.score_rules = score_rules
        self.initialize_substitution_matrix()

    def initialize_substitution_matrix(self):
        self.aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix)

    def _align_and_score(self, seq1, seq2):
        alignments = self.aligner.align(seq1[::-1], seq2[::-1])
        align_read1, align_read2 = alignments[0][0][::-1], alignments[0][1][::-1]
        scores = [self.score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
        max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
        return align_read1, align_read2, max_index

    @abstractmethod
    def align(self, seq1, seq2):
        pass

class AlignmentTransformedDefault(BaseAlignment):
    def align(self, seq1, seq2):
        align_read1, align_read2, max_index = self._align_and_score(seq1, seq2)
        return align_read1[:max_index], align_read2[:max_index]

class AlignmentTransformedLong(BaseAlignment):
    def align(self, seq1, seq2):
        align_read1, align_read2, max_index = self._align_and_score(seq1, seq2)
        align_read1 = align_read1[:max_index]
        seq1_long = seq1[len(align_read1) - align_read1.count("-"):]
        return align_read1 + seq1_long, align_read2[:max_index]

class AlignmentFactory:
    @staticmethod
    def create_alignment(transformed="CT", long_mode=False, mode='global', open_gap_score=-3, extend_gap_score=-2):
        # 检查输入有效性
        if transformed not in ["CT", "CTGA"]:
            raise ValueError(f"Invalid transformed type: {transformed}. Must be 'CT' or 'CTGA'.")

        # 根据是否为 Long 模式选择类
        if transformed == "CT":
            if long_mode:
                return AlignmentTransformedLong(
                    rule_matrix=AlignmentConfig.get_rule_matrix_ct(),
                    score_rules=AlignmentConfig.get_score_rules_ct(),
                    mode=mode,
                    open_gap_score=open_gap_score,
                    extend_gap_score=extend_gap_score
                )
            else:
                return AlignmentTransformedDefault(
                    rule_matrix=AlignmentConfig.get_rule_matrix_ct(),
                    score_rules=AlignmentConfig.get_score_rules_ct(),
                    mode=mode,
                    open_gap_score=open_gap_score,
                    extend_gap_score=extend_gap_score
                )
        elif transformed == "CTGA":
            if long_mode:
                return AlignmentTransformedLong(
                    rule_matrix=AlignmentConfig.get_rule_matrix_ctga(),
                    score_rules=AlignmentConfig.get_score_rules_ctga(),
                    mode=mode,
                    open_gap_score=open_gap_score,
                    extend_gap_score=extend_gap_score
                )
            else:
                return AlignmentTransformedDefault(
                    rule_matrix=AlignmentConfig.get_rule_matrix_ctga(),
                    score_rules=AlignmentConfig.get_score_rules_ctga(),
                    mode=mode,
                    open_gap_score=open_gap_score,
                    extend_gap_score=extend_gap_score
                )








class CorrectConfig:
    """
    配置类，用于管理严格/宽松和 CT/CTGA 的矩阵定义。
    """
    @staticmethod
    def get_mapping(strict_mode=False, transformed="CTGA"):
        """
        根据 strict_mode 和 transformed 参数返回对应的映射矩阵。
        :param strict_mode: bool, 是否严格模式
        :param transformed: str, "CT" 或 "CTGA"
        :return: dict, 映射矩阵
        """
        if transformed == "CT":
            return CorrectConfig.strict_ct_mapping() if strict_mode else CorrectConfig.relaxed_ct_mapping()
        elif transformed == "CTGA":
            return CorrectConfig.strict_ctga_mapping() if strict_mode else CorrectConfig.relaxed_ctga_mapping()
        else:
            raise ValueError(f"Invalid transformed type: {transformed}. Must be 'CT' or 'CTGA'.")

    @staticmethod
    def strict_ct_mapping():
        return {
            ('A', 'A'): 'A', ('A', 'C'): 'N', ('A', 'G'): 'N', ('A', 'T'): 'N', ('A', 'N'): 'N', ('A', '-'): 'N',
            ('C', 'A'): 'N', ('C', 'C'): 'C', ('C', 'G'): 'N', ('C', 'T'): 'N', ('C', 'N'): 'N', ('C', '-'): 'N',
            ('G', 'A'): 'N', ('G', 'C'): 'N', ('G', 'G'): 'G', ('G', 'T'): 'N', ('G', 'N'): 'N', ('G', '-'): 'N',
            ('T', 'A'): 'N', ('T', 'C'): 'C', ('T', 'G'): 'N', ('T', 'T'): 'T', ('T', 'N'): 'N', ('T', '-'): 'N',
            ('N', 'A'): 'N', ('N', 'C'): 'N', ('N', 'G'): 'N', ('N', 'T'): 'N', ('N', 'N'): 'N', ('N', '-'): 'N',
            ('-', 'A'): 'N', ('-', 'C'): 'N', ('-', 'G'): 'N', ('-', 'T'): 'N', ('-', 'N'): 'N', ('-', '-'): 'N'
        }

    @staticmethod
    def strict_ctga_mapping():
        return {
            ('A', 'A'): 'A', ('A', 'C'): 'N', ('A', 'G'): 'N', ('A', 'T'): 'N', ('A', 'N'): 'N', ('A', '-'): 'N',
            ('C', 'A'): 'N', ('C', 'C'): 'C', ('C', 'G'): 'N', ('C', 'T'): 'N', ('C', 'N'): 'N', ('C', '-'): 'N',
            ('G', 'A'): 'G', ('G', 'C'): 'N', ('G', 'G'): 'G', ('G', 'T'): 'N', ('G', 'N'): 'N', ('G', '-'): 'N',
            ('T', 'A'): 'N', ('T', 'C'): 'C', ('T', 'G'): 'N', ('T', 'T'): 'T', ('T', 'N'): 'N', ('T', '-'): 'N',
            ('N', 'A'): 'N', ('N', 'C'): 'N', ('N', 'G'): 'N', ('N', 'T'): 'N', ('N', 'N'): 'N', ('N', '-'): 'N',
            ('-', 'A'): 'N', ('-', 'C'): 'N', ('-', 'G'): 'N', ('-', 'T'): 'N', ('-', 'N'): 'N', ('-', '-'): 'N'
        }

    @staticmethod
    def relaxed_ct_mapping():
        return {('A', 'A'): 'A', ('G', 'G'): 'G', ('C', 'C'): 'C', ('T', 'T'): 'T',
                ('A', '-'): 'A', ('G', '-'): 'G', ('C', '-'): 'C', ('T', '-'): 'T',
                ('-', 'A'): 'A', ('-', 'G'): 'G', ('-', 'C'): 'C', ('-', 'T'): 'T',
                ('T', 'C'): 'C'}

    @staticmethod
    def relaxed_ctga_mapping():
        return {('A', 'A'): 'A', ('G', 'G'): 'G', ('C', 'C'): 'C', ('T', 'T'): 'T',
                ('A', '-'): 'A', ('G', '-'): 'G', ('C', '-'): 'C', ('T', '-'): 'T',
                ('-', 'A'): 'A', ('-', 'G'): 'G', ('-', 'C'): 'C', ('-', 'T'): 'T',
                ('T', 'C'): 'C', ('G', 'A'): 'G'}

class CorrectStrategy(ABC):
    def __init__(self, mapping ):
        self.mapping = mapping

    @staticmethod
    def correct_quality(str1, str2):
        str2_iter = iter(str2)
        return ''.join(',' if char == '-' else next(str2_iter, '') for char in str1)

    def mapping_correct_read(self, read1, read2, qual1, qual2):
        correct_qual1 = self.correct_quality(read1, qual1)
        correct_qual2 = self.correct_quality(read2, qual2)
        correct_qual = ''.join(max(a, b) for a, b in zip(correct_qual1, correct_qual2))
        result = [
            self.mapping[(c1, c2)] if (c1, c2) in self.mapping else (c1 if q1 >= q2 else c2)
            for c1, c2, q1, q2 in zip(read1, read2, correct_qual1, correct_qual2)
        ]
        correct_read = ''.join(result)
        correct_read1 = read1.replace("-", "N")
        correct_read2 = read2.replace("-", "N")
        return correct_read, correct_read1, correct_read2, correct_qual, correct_qual1, correct_qual2

    @abstractmethod
    def correct_read(self, read1, read2, qual1, qual2):
        pass

class DefaultCorrectStrategy(CorrectStrategy):
    def correct_read(self, read1, read2, qual1, qual2):
        return self.mapping_correct_read(read1, read2, qual1, qual2)

class LongestMatchCorrectStrategy(CorrectStrategy):
    def correct_read(self, read1, read2, qual1, qual2):
        remain_read1 = read1[len(read2):]
        read1 = read1[:len(read2)]
        correct_read, correct_read1, correct_read2, correct_qual, correct_qual1, correct_qual2 = self.mapping_correct_read(read1, read2, qual1, qual2)
        qual1_part = qual1[-len(remain_read1):] if remain_read1 else ""
        return correct_read + remain_read1, correct_read1, correct_read2, correct_qual + qual1_part, correct_qual1, correct_qual2

class CorrectStrategyFactory:
    @staticmethod
    def get_strategy(strategy_type, strict_mode=False, transformed="CTGA"):
        # 从配置类获取对应矩阵
        mapping = CorrectConfig.get_mapping(strict_mode=strict_mode, transformed=transformed)

        # 根据 strategy_type 返回对应的策略实例
        if strategy_type == "default":
            return DefaultCorrectStrategy(mapping)
        elif strategy_type == "long":
            return LongestMatchCorrectStrategy(mapping)
        else:
            raise ValueError(f"Unknown strategy type: {strategy_type}")





class PairedRead:
    def __init__(self, name, read1, read2, qual1, qual2):
        self.name = name
        self.read1 = read1
        self.read2 = read2
        self.qual1 = qual1
        self.qual2 = qual2
        self.trimmed_read1 = None
        self.trimmed_read2 = None
        self.corrected_read1 = None
        self.corrected_read2 = None
        self.corrected_read = None
        self.corrected_qual1 = None
        self.corrected_qual2 = None
        self.corrected_qual = None



class DataProcessor:
    def __init__(self, trim_strategy_type,alignment_rule, correct_strategy_type, strict_mode, transformed, min_length, processes):
        """
        初始化数据处理器。
        """

        self.trim_strategy = TrimStrategyFactory.get_strategy(trim_strategy_type)  # 用户选择修剪策略
        self.alignment = AlignmentFactory.create_alignment(alignment_rule)  # 用户选择比对规则
        self.correct_strategy = CorrectStrategyFactory.get_strategy(
            correct_strategy_type, strict_mode=strict_mode, transformed=transformed
        )  # 用户选择纠正策略
        self.min_length = min_length  # 最小长度限制
        self.processes = processes  # 并行进程数
        self.total_records = 0  # 总记录数
        self.processed_records = 0  # 已处理记录数

    def process_chunk(self, chunks):
        result = []
        for paired_read in chunks:
            # 修剪
            paired_read.trimmed_read1, paired_read.trimmed_read2 = self.trim_strategy.trim(
                paired_read.read1, paired_read.read2
            )
            # 比对
            paired_read.trimmed_read1, paired_read.trimmed_read2 = self.alignment.align(
                paired_read.trimmed_read1, paired_read.trimmed_read2
            )
            # 长度过滤
            if len(paired_read.trimmed_read1) < self.min_length or len(paired_read.trimmed_read2) < self.min_length:
                continue
            # 纠正
            paired_read.corrected_read, paired_read.corrected_read1, paired_read.corrected_read2, \
                paired_read.corrected_qual, paired_read.corrected_qual1, paired_read.corrected_qual2 = \
                self.correct_strategy.correct_read(
                    paired_read.trimmed_read1, paired_read.trimmed_read2, paired_read.qual1, paired_read.qual2
                )
            result.append(paired_read)
        return result


    def process_buffer(self, buffer):
        # 分块处理
        chunk_size = max(1, len(buffer) // self.processes)
        chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
        with Pool(processes=self.processes) as pool:
            results = pool.map(self.process_chunk, chunks)

        # 展平结果并返回给主进程
        flattened_results = [record for result in results for record in result]
        self.processed_records += len(flattened_results)
        return flattened_results

    def get_summary(self):
        """
        返回处理结果的统计信息。
        """
        return {
            "total_records": self.total_records,
            "processed_records": self.processed_records,
            "resolved_ratio": f"{(self.processed_records / self.total_records) * 100:.2f}%" if self.total_records else "N/A"
        }



class DNAProcessor:
    def __init__(self, input_file1, input_file2, output_prefix, alignment, trim_strategy, correct_strategy, strict_mode, transformed, chunk_size=100000, min_length=50, processes=24):
        self.file_reader = FileReader(input_file1, input_file2)  # 文件读取器
        self.output_prefix = output_prefix
        self.data_processor = DataProcessor(
            alignment_rule=alignment,
            trim_strategy_type=trim_strategy,
            correct_strategy_type=correct_strategy,
            strict_mode=strict_mode,
            transformed=transformed,
            min_length=min_length,
            processes=processes,
        )
        self.chunk_size = chunk_size  # 缓冲区大小
        self.buffer = []  # 初始化缓冲区

    def write_sequence(self, paired_read, output_handles):
        output_handles[0].write(f"{paired_read.name}\n")
        output_handles[0].write(f"{paired_read.corrected_read}\n")
        output_handles[0].write("+\n")
        output_handles[0].write(f"{paired_read.corrected_qual}\n")
        output_handles[1].write(f"{paired_read.name}\n")
        output_handles[1].write(f"{paired_read.corrected_read1}\n")
        output_handles[1].write("+\n")
        output_handles[1].write(f"{paired_read.corrected_qual1}\n")
        output_handles[2].write(f"{paired_read.name}\n")
        output_handles[2].write(f"{paired_read.corrected_read2}\n")
        output_handles[2].write("+\n")
        output_handles[2].write(f"{paired_read.corrected_qual2}\n")

    def process(self):
        """
        主处理逻辑。
        """
        output_handles = [
            open(self.output_prefix + ".fq", "w"),
            open(self.output_prefix + "_cut_f1.fq", "w"),
            open(self.output_prefix + "_cut_f2.fq", "w"),
        ]

        try:
            for record in self.file_reader:
                self.buffer.append(record)
                self.data_processor.total_records += 1

                if len(self.buffer) >= self.chunk_size:
                    results = self.data_processor.process_buffer(self.buffer)
                    for paired_read in results:
                        self.write_sequence(paired_read, output_handles)
                    self.buffer = []

            if self.buffer:
                results = self.data_processor.process_buffer(self.buffer)
                for paired_read in results:
                    self.write_sequence(paired_read, output_handles)

            # 打印统计结果
            summary = self.data_processor.get_summary()
            print(
                f"Total Records: {summary['total_records']}, "
                f"Processed Records: {summary['processed_records']}, "
                f"Resolved Ratio: {summary['resolved_ratio']}"
            )
        finally:
            # 关闭文件资源
            self.file_reader.close()
            for handle in output_handles:
                handle.close()




def main():
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="DNA Sequencing Processor Tool")
    parser.add_argument("-i1", "--fq1", required=True, help="Path to input FASTQ file 1")
    parser.add_argument("-i2", "--fq2", required=True, help="Path to input FASTQ file 2")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("-a", "--alignment", choices=["CT", "CTGA"], default="CTGA", help="Alignment rule type")
    parser.add_argument("-t", "--trim_strategy", choices=["default", "advanced"], default="default", help="Trimming strategy type")
    parser.add_argument("-c", "--correct_strategy", choices=["default", "long"], default="default", help="Correction strategy type")
    parser.add_argument("-s", "--strict_mode", action="store_true", help="Enable strict mode for correction")
    parser.add_argument("-f", "--transformed", choices=["CT", "CTGA"], default="CTGA", help="Transformation type for correction")
    parser.add_argument("-m", "--min", type=int, default=50, help="Minimum length for processed reads")
    parser.add_argument("-p", "--processes", type=int, default=20, help="Number of processes for parallel processing")
    parser.add_argument("-b", "--chunk_size", type=int, default=100000, help="Buffer size for processing chunks")

    # 解析参数
    args = parser.parse_args()

    # 输出参数信息以供调试
    print("Running DNA Processor with the following parameters:")
    print(f"Input File 1: {args.fq1}")
    print(f"Input File 2: {args.fq2}")
    print(f"Output Prefix: {args.output_prefix}")
    print(f"Alignment Rule: {args.alignment}")
    print(f"Trimming Strategy: {args.trim_strategy}")
    print(f"Correction Strategy: {args.correct_strategy}")
    print(f"Strict Mode: {'Enabled' if args.strict_mode else 'Disabled'}")
    print(f"Transformation Type: {args.transformed}")
    print(f"Minimum Length: {args.min}")
    print(f"Number of Processes: {args.processes}")
    print(f"Chunk Size: {args.chunk_size}")

    # 初始化 DNAProcessor 并执行
    processor = DNAProcessor(
        input_file1=args.fq1,
        input_file2=args.fq2,
        output_prefix=args.output_prefix,
        alignment=args.alignment,
        trim_strategy=args.trim_strategy,
        correct_strategy=args.correct_strategy,
        strict_mode=args.strict_mode,
        transformed=args.transformed,
        chunk_size=args.chunk_size,
        min_length=args.min,
        processes=args.processes,
    )

    # 执行处理逻辑
    processor.process()

if __name__ == "__main__":
    main()




