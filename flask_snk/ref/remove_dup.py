import pysam
import Levenshtein
import argparse
import gzip
from collections import defaultdict




def get_pos_read(input_bam):
    chr_pos_trend=defaultdict(list)
    with pysam.AlignmentFile(input_bam, "rb") as input_bamfile:
        for read in input_bamfile:
            if read.is_unmapped:
                continue
            chrom = read.reference_name
            pos1 = read.reference_start
            pos2 = read.reference_end
            is_reverse = read.is_reverse
            if is_reverse:
                chr_pos_trend[chrom + "_" + str(pos2) + "_" + str(int(is_reverse))].append(read.qname)
            else:
                chr_pos_trend[chrom + "_" + str(pos1) + "_" + str(int(is_reverse))].append(read.qname)

    return chr_pos_trend

def get_pair_fq(fq1,fq2):
    fq_quality={}
    with gzip.open(fq1,"rb")as file1,gzip.open(fq2,"rb")as file2:
        for i,line in enumerate(zip(file1,file2)):
            t1 = line[0].decode('utf-8')
            t2 = line[1].decode('utf-8')
            if i %4 ==0:
                name=t1.strip().split()[0]
            if i %4==1:
                seq1=t1.strip()
                seq2=t2.strip()
            if i %4==3:
                q1=t1.strip()
                q2=t2.strip()
                q_score = int ((sum(ord(char) for char in q1) + sum(ord(char) for char in q2) ) / (len(q1)+len(q2)) * 10000)
                fq_quality[name[1:]]=(seq1,seq2,q_score)

    return fq_quality

def final_dup_from_list(sorted_list,temp_dir,dis):
    result_list=[sorted_list[0]]
    for read in sorted_list:
        temp_seq=temp_dir[read]
        temp_seq1=temp_seq[0]
        temp_seq2=temp_seq[1]

        temp_re = read
        for re_name in result_list:
            re_seq=temp_dir[re_name]
            re_seq1=re_seq[0]
            re_seq2=re_seq[1]
            if Levenshtein.distance(temp_seq1,re_seq1)<dis and Levenshtein.distance(temp_seq2,re_seq2)<dis:
                temp_re=""
                break
        if temp_re != "":
            result_list.append(temp_re)


    return result_list

def write_dup_bam(chr_pos_trend,fq_quality,input_bam,output_bam,dis):
    dup_list_name=[]
    for name_list in chr_pos_trend.values():
        temp_list=[(name,fq_quality[name][-1]) for name in name_list]
        sorted_list = [temp[0] for temp in sorted(temp_list, key=lambda x: x[-1], reverse=True)]
        # 优化一下，提前获取小字典
        temp_dir=  {key: fq_quality[key] for key in name_list}
        result_list=final_dup_from_list(sorted_list,temp_dir,dis)
        dup_list_name.extend(result_list)

    dup_list_name=set(dup_list_name)
    with pysam.AlignmentFile(input_bam, "rb") as input_bamfile, pysam.AlignmentFile(output_bam, "wb",
                                                                                    header=input_bamfile.header) as output_bamfile:
        for read in input_bamfile:
            if read.qname in dup_list_name:
                output_bamfile.write(read)




def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("--input", type=str, help="Input fq1 file")
    parser.add_argument("--output", type=str, help="Output w file")
    parser.add_argument("--fq1", type=str, help="Input fq1 file")
    parser.add_argument("--fq2", type=str, help="Input fq2 file")

    parser.add_argument("--dis", type=int, default=5, help="default: 5")

    args = parser.parse_args()

    chr_pos_trend = get_pos_read(args.input)
    fq_quality = get_pair_fq(args.fq1,args.fq2)
    write_dup_bam(chr_pos_trend, fq_quality, args.input,args.output,args.dis)


if __name__ == '__main__':
    main()


