#! /usr/bin/env python3
import argparse


def determine_homopolymer(ref, output, homopolymer_length):
    chromosome = base = start = end = -1
    with open(ref, 'r') as input_file:
        with open(output, 'w') as write_file:
            for line in input_file:
                if '>' in line:
                    if end - start > homopolymer_length:
                        write_file.write('\t'.join(map(str, [chromosome, start, end, end - start])))
                    chromosome = line.strip().split()[0][1:]
                    base = -1
                    start = end = 0
                else:
                    for b in line.strip():
                        if b != base:
                            if end - start > homopolymer_length:
                                write_file.write('\t'.join(map(str, [chromosome, start, end, end - start])) + "\n")
                            start = end
                            base = b
                        end += 1

            if end - start > homopolymer_length:
                write_file.write('\t'.join(map(str, [chromosome, start, end, end - start])) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', help='full path to reference fasta')
    parser.add_argument('output_file', help='name of outputfile')
    parser.add_argument(
        '--homopolymer_len',
        default=3,
        type=int,
        help='minimum length of homopolymer to consider in output [default 3]'
    )
    args = parser.parse_args()
    determine_homopolymer(args.reference, args.output_file, args.homopolymer_len)
