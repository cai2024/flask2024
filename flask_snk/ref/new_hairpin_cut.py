import Levenshtein
import gzip
import math
from functools import partial
import multiprocessing
from Bio.Align.substitution_matrices import Array
from Bio.Align import PairwiseAligner
import numpy as np
import argparse


class Alignment:
    alphabet = "ACGTN"
    rule_matrix_1 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    rule_matrix_2 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    def __init__(self, sub_matrix, mode='global', open_gap_score=-3, extend_gap_score=-2):
        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score
        if sub_matrix == 'rule_matrix_1':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_1)
        elif sub_matrix == 'rule_matrix_2':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_2)
        else:
            raise ValueError(f"Invalid sub_matrix argument: {sub_matrix}")
        self.aligner = aligner
        self.type = sub_matrix

    def align(self, seq1, seq2):
        alignments = self.aligner.align(seq1[::-1], seq2[::-1])
        align_read1, align_read2 = alignments[0][0][::-1], alignments[0][1][::-1]
        if self.type == 'rule_matrix_1':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): -20, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
            }
            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]

        elif self.type == 'rule_matrix_2':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): 1, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
            }
            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]
        else:
            raise ValueError(f"Invalid sub_matrix argument: {self.type}")

        return  align_read1, align_read2


class DNAProcessor:
    def __init__(self, sub_matrix):
        if sub_matrix == 'rule_matrix_1':
            self.correct_seq = self._correct_seq_option1
        elif sub_matrix == 'rule_matrix_2':
            self.correct_seq = self._correct_seq_option2
        else:
            raise ValueError(f"Invalid sub_matrix argument: {sub_matrix}")

    def _correct_seq_option1(self,seq1, seq2, qual1, qual2):
        correct_qual1 = []
        correct_qual2 = []
        correct_qual = []
        correct_read1 = []
        correct_read2 = []
        correct_read = []
        temp_r1 = 0
        temp_r2 = 0
        for i, (char1, char2) in enumerate(zip(seq1, seq2)):
            q1 = qual1[i + temp_r1]
            q2 = qual2[i + temp_r2]
            if char1 == char2:
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(max(q1, q2))
            elif char1 == "-":
                temp_r1 -= 1
                correct_read1.append("N")
                correct_read2.append(char2)
                correct_read.append(char2)
                correct_qual1.append(",")
                correct_qual2.append(q2)
                correct_qual.append(q2)
            elif char2 == "-":
                temp_r2 -= 1
                correct_read1.append(char1)
                correct_read2.append("N")
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(",")
                correct_qual.append(q1)

            elif char1=="T" and char2=="C":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append("C")
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(max(q1, q2))

            elif q1 > q2 and seq1 != "T":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(q1)
            elif q2 > q1 :
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char2)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(q2)
            else:
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append("N")
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(",")
        return "".join(correct_read1), "".join(correct_read2), "".join(correct_read), "".join(correct_qual1), "".join(
            correct_qual2), "".join(correct_qual)

    def _correct_seq_option2(self,seq1, seq2, qual1, qual2):
        correct_qual1 = []
        correct_qual2 = []
        correct_qual = []
        correct_read1 = []
        correct_read2 = []
        correct_read = []
        temp_r1 = 0
        temp_r2 = 0
        for i, (char1, char2) in enumerate(zip(seq1, seq2)):
            q1 = qual1[i + temp_r1]
            q2 = qual2[i + temp_r2]
            if char1 == char2:
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(max(q1, q2))
            elif char1 == "-":
                temp_r1 -= 1
                correct_read1.append("N")
                correct_read2.append(char2)
                correct_read.append(char2)
                correct_qual1.append(",")
                correct_qual2.append(q2)
                correct_qual.append(q2)
            elif char2 == "-":
                temp_r2 -= 1
                correct_read1.append(char1)
                correct_read2.append("N")
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(",")
                correct_qual.append(q1)

            elif char1=="T" and char2=="C":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append("C")
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(max(q1, q2))
            elif char1=="G" and char2=="A":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append("G")
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(max(q1, q2))


            elif q1 > q2 and seq1 != "T":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char1)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(q1)
            elif q2 > q1 and seq2 != "A":
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append(char2)
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(q2)
            else:
                correct_read1.append(char1)
                correct_read2.append(char2)
                correct_read.append("N")
                correct_qual1.append(q1)
                correct_qual2.append(q2)
                correct_qual.append(",")
        return "".join(correct_read1), "".join(correct_read2), "".join(correct_read), "".join(correct_qual1), "".join(
            correct_qual2), "".join(correct_qual)


    def open_input_file(self, input_file):
        if input_file.endswith('.gz'):
            in_handle = gzip.open(input_file, "rt")
        else:
            in_handle = open(input_file)
        return in_handle
    def write_sequence(self, seq, file_handle, index):
        file_handle.write(f"{seq[-1]}\n")
        file_handle.write(f"{seq[index]}\n")
        file_handle.write("+\n")
        file_handle.write(f"{seq[index + 3]}\n")
    def complement_dna(self,dna):
        complement = str.maketrans('ATCG', 'TAGC')
        reverse = dna[::-1].translate(complement)
        return reverse
    def check(self,str1, str2):
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
    def trim_overlap(self,read1, read2):
        a1 = self.check(read1, self.complement_dna(read2))
        if a1 > 5:
            return read1[:-int(a1 / 2)], read2[:-int(a1 / 2)]
        else:
            a2 = self.check(self.complement_dna(read2), read1)
            if a2 > min(len(read2), len(read1)) * 0.8 and a2 > 5:
                return read1[:-int(a2 / 2)], read2[:-int(a2 / 2)]
        return read1, read2

    def process_chunk(self, chunks, alignment, min_lenth=50):
        result = []
        for record in chunks:
            read1, read2 = self.trim_overlap(record[0], record[1])
            seq1, seq2 = alignment.align(read1, read2)
            if len(seq1) < min_lenth:
                continue
            f1,f2,f3,f4,f5,f6 = self.correct_seq(seq1, seq2, record[2], record[3])
            result.append([f1,f2,f3,f4,f5,f6,record[-1]])
        return result

    def process_records(self, buffer, alignment, min_lenth=50, processes=24):
        pool = multiprocessing.Pool(processes=processes)
        chunk_size = math.ceil(len(buffer) / processes)
        chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
        partial_process_chunk = partial(self.process_chunk, alignment=alignment, min_lenth=min_lenth)
        results = pool.map(partial_process_chunk, chunks)
        pool.close()
        pool.join()
        return [record for result in results for record in result]


def cache_and_process(input_file1, input_file2, output_file, alignment,dnaprocessor, min_lenth=50, processes=24, chunk_size=1000000):
    cut_f1_path = output_file +'_cut_f1.fq'
    cut_f2_path = output_file +'_cut_f2.fq'
    in_handle1 = dnaprocessor.open_input_file(input_file1)
    in_handle2 = dnaprocessor.open_input_file(input_file2)

    report_n1 = 0
    with (in_handle1, in_handle2, open(output_file+".fq", "w") as out_handle, open(cut_f1_path, "w") as f1, open(cut_f2_path, "w") as f2):
        buffer = []
        for line_num, (line1, line2) in enumerate(zip(in_handle1, in_handle2)):
            if line_num % 4 == 0:
                record_name = line1.strip().split()[0]
            elif line_num % 4 == 1:
                record_seq1 = line1.strip()
                record_seq2 = line2.strip()
            elif line_num % 4 == 3:
                record_qual1 = line1.strip()
                record_qual2 = line2.strip()
                buffer.append([record_seq1,record_seq2,record_qual1,record_qual2,record_name])
            if len(buffer) == chunk_size:
                print("process reads : " + str(int((line_num + 1) / 4)))
                temp = dnaprocessor.process_records(buffer, alignment, min_lenth, processes)
                report_n1 = report_n1 + len(temp)
                for seq in temp:
                    dnaprocessor.write_sequence(seq, f1, 0)
                    dnaprocessor.write_sequence(seq, f2, 1)
                    dnaprocessor.write_sequence(seq, out_handle, 2)
                buffer = []
        if buffer:
            temp = dnaprocessor.process_records(buffer, alignment, min_lenth, processes)
            report_n1 = report_n1 + len(temp)
            for seq in temp:
                dnaprocessor.write_sequence(seq, f1, 0)
                dnaprocessor.write_sequence(seq, f2, 1)
                dnaprocessor.write_sequence(seq, out_handle, 2)


    message = "fq1: " + input_file1 + "\n" + "fq2: " + input_file2 + "\n" + "min_read_length: " + str(
        min_lenth) + "\n" + "input_reads: " + str(int((line_num + 1) / 4)) + "\n" + "resolved_reads: " + str(
        report_n1) + "\n" + "resolved ratio: " + str("{:.2f}%".format(report_n1 / int((line_num + 1) / 4) * 100))
    print(message)
    return message





def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("--fq1", type=str, help="Input fq1 file")
    parser.add_argument("--fq2", type=str, help="Input fq2 file")
    parser.add_argument("--outfile", type=str,default="restored_seq", help="Output w file")
    parser.add_argument("--rule", type=int, default=2, help="default: 2")
    parser.add_argument("--min", type=int, default=50, help="Position parameter")
    parser.add_argument("--parallel", type=int, default=24, help="default: 24")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="default: 1000000")

    args = parser.parse_args()

    if args.rule == 2:
        alignment = Alignment("rule_matrix_2")
        dnaprocessor=DNAProcessor("rule_matrix_2")
        mylog = cache_and_process(args.fq1, args.fq2, args.outfile, alignment,dnaprocessor, args.min, args.parallel, args.chunk_size)
    else:
        alignment = Alignment("rule_matrix_1")
        dnaprocessor = DNAProcessor("rule_matrix_1")
        mylog = cache_and_process(args.fq1, args.fq2, args.outfile, alignment, dnaprocessor, args.min, args.parallel,
                                  args.chunk_size)



    print(mylog)


if __name__ == '__main__':
    main()









