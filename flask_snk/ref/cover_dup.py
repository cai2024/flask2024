import argparse


def get_log(cov, dep, bed_dir1, bed_dirn, chr_dir, spikein, output):
    with open(output, "a") as file:
        c1 = 0
        cn = 0
        for key, value in cov.items():
            if key in spikein:
                print(key + "标品cover:" + str(value), file=file)
            else:
                c1 += bed_dir1[key]
                cn += chr_dir[key]
        print("基因组cover：" + str(c1 / cn), file=file)

        c1 = 0
        cn = 0
        for key, value in dep.items():
            if key in spikein:
                print(key + "标品depth:" + str(value), file=file)
            else:
                c1 += bed_dir1[key]
                cn += bed_dirn[key]
        print("基因组dep：" + str(cn / c1), file=file)

        print("", file=file)
        print("", file=file)
        print("chr的cover", file=file)
        for i, j in cov.items():
            print(i + " : " + str(j), file=file)
        print("", file=file)
        print("", file=file)
        print("chr的depth", file=file)
        for i, j in dep.items():
            print(i + " : " + str(j), file=file)








if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fai",
        default="/data/wangzc/cgs/new_bis/new_file.fa.fai",
        help="fai file",
        type=str,
    )

    parser.add_argument(
        "--bed",
        help="bed file",
        type=str,
    )

    parser.add_argument(
        "--spikein",
        default="lambda,mc,G5hmc,puc19",
        help="barcodes file",
        type=str,
    )

    parser.add_argument(
        "--output",
        help="output file",
        default="N",
        type=str,
    )
    args = parser.parse_args()

    if args.output=="N":
        args.output=args.bed+".cov.log"

    spikein = args.spikein.split(',')

    chr_dir = {}
    with open(args.fai) as file:
        for line in file:
            line = line.split()
            chr_dir[line[0]] = int(line[1])

    bed_dir1 = {key: 0 for key in chr_dir.keys()}
    bed_dirn = {key: 0 for key in chr_dir.keys()}
    cov = {}
    dep = {}

    with open(args.bed) as file:
        for line in file:
            line = line.split()
            line[3]=int(float(line[3]))
            if int(line[3]) == 0:
                continue
            elif int(line[3]) == 1:
                bed_dir1[line[0]] += (int(line[2]) - int(line[1]))
                bed_dirn[line[0]] += (int(line[2]) - int(line[1]))
            else:
                bed_dir1[line[0]] += (int(line[2]) - int(line[1]))
                bed_dirn[line[0]] += ((int(line[2]) - int(line[1])) * int(line[3]))

    for i, j in bed_dir1.items():
        if chr_dir[i] != 0:
            cov[i] = j / chr_dir[i]
        else:
            cov[i] = 0
    for i, j in bed_dirn.items():
        if bed_dir1[i] != 0:
            dep[i] = j / bed_dir1[i]
        else:
            dep[i] = 0

    get_log(cov,dep,bed_dir1,bed_dirn,chr_dir,spikein,args.output)




