import argparse
import subprocess



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        default="output.log",
        type=str,
    )
    parser.add_argument(
        "--hmc",
        help="bed file",
        type=str,
        default="hmc.bed",
    )
    parser.add_argument(
        "--mc",
        help="bed file",
        type=str,
        default="mc.bed",
    )

    args = parser.parse_args()

    if args.output=="output.log":
        args.output=args.bed+".cov.log"
    spikein = args.spikein.split(',')
    ref={"lambda":48502,
         "mc":897,
         "G5hmc":335,
         "puc19":2686,
    }

    cover_base = {key: 0 for key in spikein}
    depth_base = {key: 0 for key in spikein}

    with open(args.bed) as file:
        for line in file:
            line = line.split()
            if int(float(line[3])) != 0 and line[0] in spikein:
                cover_base[line[0]] += (int(line[2]) - int(line[1]))
                depth_base[line[0]] += ((int(line[2]) - int(line[1])) * int(float(line[3])))
    result=[]

    if args.hmc=="hmc.bed":
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
                command1 = f"""awk '$7 == "CpG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.hmc}"""
                command2 = f"""awk '$7 == "CHG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.hmc}"""
                command3 = f"""awk '$7 == "CHH" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.hmc}"""
                temp1 = subprocess.run(command1, shell=True, text=True, capture_output=True)
                temp2 = subprocess.run(command2, shell=True, text=True, capture_output=True)
                temp3 = subprocess.run(command3, shell=True, text=True, capture_output=True)
                cpg = 0 if temp1.stdout.strip()=="" else round(float(temp1.stdout.strip()),4)
                chg = 0 if temp2.stdout.strip() == "" else round(float(temp2.stdout.strip()), 4)
                chh = 0 if temp3.stdout.strip() == "" else round(float(temp3.stdout.strip()), 4)
                command1 = f"""awk '$7 == "CpG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.mc}"""
                command2 = f"""awk '$7 == "CHG" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.mc}"""
                command3 = f"""awk '$7 == "CHH" && $1 == "{key}" {{ sum += $4; count++ }} END {{ if (count > 0) print sum / count }}' {args.mc}"""
                temp1 = subprocess.run(command1, shell=True, text=True, capture_output=True)
                temp2 = subprocess.run(command2, shell=True, text=True, capture_output=True)
                temp3 = subprocess.run(command3, shell=True, text=True, capture_output=True)
                mc_cpg = 0 if temp1.stdout.strip()=="" else round(float(temp1.stdout.strip()),4)
                mc_chg = 0 if temp2.stdout.strip() == "" else round(float(temp2.stdout.strip()), 4)
                mc_chh = 0 if temp3.stdout.strip() == "" else round(float(temp3.stdout.strip()), 4)


                result.append(f'{cover}|{depth}\t{cpg}|{chg}|{chh}\t{mc_cpg}|{mc_chg}|{mc_chh}')
        print('\t'.join(result))









