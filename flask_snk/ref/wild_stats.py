import argparse

def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")
    parser.add_argument("--fq1", type=str, help="Input fq1 file")
    parser.add_argument("--fq2", type=str, help="Input fq2 file")
    parser.add_argument("--output", type=str, help="default: 2")
    args = parser.parse_args()
    if not args.output:
        args.output = args.fq1 + ".stats"
    my_dir={}

    with open(args.fq1, 'r') as file1, open(args.fq2, 'r') as file2, open(args.output,"w")as file3:
        for i,(line1, line2) in enumerate(zip(file1, file2)):
            line1 = line1.strip()
            line2 = line2.strip()
            if i %4==1:
                for char1, char2 in zip(line1, line2):
                    combination = (char1, char2)
                    if combination in my_dir:
                        my_dir[combination] += 1
                    else:
                        my_dir[combination] = 1

        sorted_items = sorted(my_dir.items(), key=lambda item: item[1], reverse=True)
        for key, value in sorted_items:
            print(f"{key}: {value}", file=file3)
        wild_c = my_dir[('T', 'C')]/(my_dir[('T', 'C')]+my_dir[('C', 'C')])
        wild_g = my_dir[('G', 'A')] / (my_dir[('G', 'A')] + my_dir[('G', 'G')])
        print("wild_c_ratio:" + str(wild_c), file=file3)
        print("wild_g_ratio:" + str(wild_g), file=file3)



if __name__ == '__main__':
    main()

