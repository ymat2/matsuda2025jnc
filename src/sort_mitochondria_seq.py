import argparse
from pathlib import Path

"""
Sort mitochondrial sequence to align AP003321 used in MitoToolPy.
D-loop region (15,554-) appears in front and there is a deletion before 859th nucleotide of d-loop.
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    mito_seq = fasta2dict(args.input)
    for h in mito_seq.keys():
        mito_seq[h] = allocate_dloop(mito_seq[h], 15554)
        mito_seq[h] = insert_deletion(mito_seq[h], 859)
    write_fasta(mito_seq, args.output)


def fasta2dict(fasta: Path) -> dict:
    seq_dict = {}
    with open(fasta) as f:
        tmp_header = ""
        for line in f:
            if line[0] == ">":
                tmp_header = line[1:]
                seq_dict[tmp_header] = ""
            else:
                seq_dict[tmp_header] += line.rstrip("\n")
    return seq_dict


def write_fasta(seq_dict: dict, fasta: Path):
    with open(fasta, "w") as f:
        for header, seq in seq_dict.items():
            f.write(">"+header+seq+"\n")


def allocate_dloop(seq: str, pos: int) -> str:
    block_front = seq[:pos-1]
    block_rear = seq[pos-1:]
    new_seq = block_rear + block_front
    return new_seq


def insert_deletion(seq: str, pos: int) -> str:
    block_front = seq[:pos-1]
    block_rear = seq[pos-1:]
    new_seq = block_front + "-" + block_rear
    return new_seq


if __name__ == "__main__":
    main()
