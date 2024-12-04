import argparse
import collections
import pysam





def call_bam_xm(bam_file,bed_file):
    update_map = {
        "Z": (0, 1, 0, 0), "z": (0, 0, 0, 1), "X": (0, 1, 0, 0), "x": (0, 0, 0, 1),
        "H": (0, 1, 0, 0), "h": (0, 0, 0, 1), "U": (0, 1, 0, 0), "u": (0, 0, 0, 1),
        "Y": (1, 0, 0, 0), "y": (0, 0, 1, 0), "V": (1, 0, 0, 0), "v": (0, 0, 1, 0),
        "W": (1, 0, 0, 0), "w": (0, 0, 1, 0), "N": (1, 0, 0, 0), "n": (0, 0, 1, 0),
    }
    bed_dict = collections.defaultdict(lambda: (0, 0, 0, 0))
    bam = pysam.AlignmentFile(bam_file, "r")
    out_file = open(bed_file , "w")

    for read in bam.fetch():
        if read.is_unmapped:
            continue
        chrom = read.reference_name
        xm_value = read.get_tag("XM")
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)

        for query_pos, ref_pos, ref_base in aligned_pairs:
            if query_pos is None:
                continue
            if xm_value[query_pos] in update_map:
                context = "CpG" if xm_value[query_pos] in "ZzYy" else "CHG" if xm_value[query_pos] in "XxVv" else "CHH" if xm_value[query_pos] in "HhWw" else "CN"
                key = (chrom, ref_pos, context)
                bed_dict[key] = tuple(map(sum, zip(bed_dict[key], update_map[xm_value[query_pos]])))

    for key, counts in bed_dict.items():
        chrom, ref_pos, context = key
        hmc_counts = counts[0]
        mc_counts = counts[1]
        unhmc_counts = counts[2]
        C_counts = counts[3]
        if (counts[0] + counts[1] + counts[3]) == 0:
            out_file.write(f"{chrom}\t{ref_pos+1}\t{ref_pos+1}\t{0:.4f}\t{0:.4f}\t{hmc_counts}\t{mc_counts}\t{unhmc_counts}\t{C_counts}\t{context}\n")
        else:
            hmc_ratio = counts[0] / (counts[0] + counts[1] + counts[3])
            mc_ratio = counts[1] / (counts[0] + counts[1] + counts[3])
            out_file.write(f"{chrom}\t{ref_pos+1}\t{ref_pos+1}\t{hmc_ratio:.4f}\t{mc_ratio:.4f}\t{hmc_counts}\t{mc_counts}\t{unhmc_counts}\t{C_counts}\t{context}\n")
    bam.close()
    out_file.close()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--parallel",
        help="Number of parallel threads to use for processing. Default: 20",
        default=20,
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
        "--bed",
        help="output file",
        type=str,
    )



    args = parser.parse_args()
    call_bam_xm(args.bam, args.bed)




if __name__ == '__main__':
    main()





