#! /usr/bin/env python3

import argparse

def split_parse_genome(reference, contig):
    with open(reference, 'r') as reference_file:
        with open(f"{contig}.fa", 'w') as write_file:
            for line in reference_file:
                if ">" in line:
                    if contig == (line.split()[0].strip(">")):
                        write_file.write(line)
                        corrrect_contig = True
                    else:
                        corrrect_contig = False
                elif corrrect_contig:
                    write_file.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', help='path to fasta of reference genome')
    parser.add_argument('contig', help='contig for output')
    args = parser.parse_args()
    split_parse_genome(args.reference, args.contig)
