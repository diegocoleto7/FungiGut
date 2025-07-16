#!/usr/bin/env python3
import sys

def parse_fasta_lengths(fasta_path):
    lengths = {}
    with open(fasta_path) as f:
        seq_id = None
        seq_len = 0
        for line in f:
            if line.startswith('>'):
                if seq_id is not None:
                    lengths[seq_id] = seq_len
                header = line[1:].strip().split()[0]
                seq_id = header
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_id is not None:
            lengths[seq_id] = seq_len
    return lengths

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(1)
    fasta = sys.argv[1]
    lengths = parse_fasta_lengths(fasta)
    for acc, L in lengths.items():
        print(f"{acc}\t{L}")
