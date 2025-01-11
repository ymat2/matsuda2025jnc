import argparse
# import gzip
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    lines, pops = parse_plink(args.input)
    write_treemix(lines, pops, args.output)


def parse_plink(plink: Path) -> (dict, list):
    rs2pop = dict()
    pops = list()

    with open(plink) as f:
        next(f)
        for line in f:
            line = line.rstrip().split()
            rs = line[1]
            pop = line[2]
            mac = line[6]
            total = line[7]

            if pop not in pops:
                pops.append(pop)

            if rs not in rs2pop:
                rs2pop[rs] = dict()
            if pop not in rs2pop[rs]:
                rs2pop[rs][pop] = str(mac)+","+str(int(total)-int(mac))

    return(rs2pop, pops)


def write_treemix(pl: dict, pops: list, tm: Path) -> None:
    with open(tm, "w") as f:
        f.write("\t".join(pops)+"\n")
        for rs in pl.values():
            rss = [rs[pop] for pop in pops]
            f.write("\t".join(rss)+"\n")


if __name__ == "__main__":
    main()
