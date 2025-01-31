#! /usr/bin/env python3
"""
Author: Joyce van der Sel
Create a altered fasta file.

python3 [script.py] [fasta.fa] [alt_file] [start-stop position] [output.fa] [log file name and sequence header] [depth_file] [sample_name_hap]

Above is the command needed to run this file.
in the alt file the strucure should be like: location_ref_alt
the location should be split by a dash (-)
It checks if alt_file has an header
"""

# import statements
import sys
import re
import os


def parse_fasta(fasta):
    """Splits fasta per chromosome and selects one (chr5)

    :param fasta: list; split per chromosome
    :return: string; one chromosoom (5 in this caase)
    """
    for chrom in fasta:
        seq = chrom.split('\n')
        if seq[0] == 'CHR5':  # here you could select a different chromosome
            seq = ''.join(seq[1:])
    print('Chrom selected', file=sys.stderr)
    return seq

def parse_depth_file(depth_file):
    """
    :param depth_file: list; per position
    :return: dictionary;
    """
    dic = {}

    for position in depth_file:
        pos = position.split("\t")
        if pos[2] not in dic:
            dic[pos[2]]=[]
            dic[pos[2]].append(int(pos[1])-1)
        else:
            dic[pos[2]].append(int(pos[1])-1)

    return dic

def replace_low_depth_nuc(dictionary, seq):
    """
    Replace nucleotide with N if depth < 3

    :param dictionary; dict; key = depth, value = POS
    :param seq:list; Every item in the list being one nucleotide of a DNA string
    """

    for key,value in dictionary.items():
        if int(key) < 3:
           for position in value:
               seq[position] = "N"
    return seq

def snp_change_fasta(seq, pos, ref, alt):
    """Changes one position in the sequence

    :param seq:list; Every item in the list being one nucleotide of a DNA string
    :param pos: int; location of the beginning of the deletion
    :param ref: string; can be one or multiple letters
    :param alt: string; can be one or multiple letters
    :return: list; altered DNA sequence

    The nucleotide of the reference allele get checked with the reference
    nucleotide at that position. Then the nucleotide item in the list gets
    replaced with a new nucleotide

    !It is important to know that the reference allele can only be of lenght 1
    otherwise it will cause an error!
    """
    if len(ref) == 1 and seq[pos] == ref:  # extra security step
        seq[pos] = alt
        print(f'{pos} has been turned from {ref} to {alt}')
    elif seq[pos] == 'N':
        print(f'depth on {pos} lower than 3')
    else:
        print(f'error in {pos}, the length of the  is {len(ref)},'
              f'the base of the reference alle is {ref} the of the fasta '
              f'nucleotide is {seq[pos]}')
    return seq


def insertion_fasta(seq, pos, ref, alt):
    """ Creates insertions in the sequence

    :param seq:list; Every item in the list being one nucleotide of a DNA string
    :param pos: int; location of the beginning of the deletion
    :param ref: sting; can be one or multiple letters
    :param alt: string; can be one or multiple letters
    :return: list; altered DNA sequence

    The nucleotide of the reference allele get checked with the reference
    nucleotide at that position. Then the nucleotide item in the list gets
    replaced with a string of the new nucleotides

    !It is important to know that the reference allele can only be of lenght 1
    otherwise it will cause an error!
    """
    if ref == seq[pos]:
        seq[pos] = alt
        print(f'{pos} has been turned from {ref} to {alt}')
    else:
        print(
            f'The reference ({ref[0]}) is different than in the fasta ({seq[pos]}) (or vise versa)')
    return seq


def deletion_fasta(seq, pos, ref, alt):
    """ Creates deletions in the sequence

    :param seq:list; Every item in the list being one nucleotide  of a DNA string
    :param pos: int; location of the beginning of the deletion
    :param ref: sting; can be one or multiple letters
    :param alt: string; can be one or multiple letters
    :return: list; altered DNA sequence

    First the nucleotides get removed then the alternative alleles get added.
    list items remain empty to maintain indexing.
    """
    ref = list(ref)
    alt = list(alt)
    if ref[0] == seq[pos]:
        for i in range(len(ref)):
            seq[pos + i] = ''
        for j in range(len(alt)):
            seq[pos + j] = alt[j]
        print(f'{pos} has been turned from {ref} to {alt}')
    else:
        print(
            f'The reference ({ref[0]}) is different than in the fasta ({seq[pos]}) (or vice versa)')
    return seq



def alter_sequence(seq, positions):
    """It goed through the positions given and uses the alterations functions

    :param seq:list; Every item in the list being one nucleotide of a DNA string
    :param posistions: list; every item an combination of location,'ref','alt'
    :return: list; altered DNA sequence
    """

    for altering in positions:
        altering = altering.split('\t')
        altering = ''.join(altering[0])
        altering = altering.split('_')
        pos = int(altering[0]) - 1
        ref = altering[1]
        alt = altering[2]


        if len(ref) == len(alt):
            seq = snp_change_fasta(seq, pos, ref, alt)
        elif len(ref) < len(alt):
            seq = insertion_fasta(seq, pos, ref, alt)
        elif len(ref) > len(alt):
            seq = deletion_fasta(seq, pos, ref, alt)
    return seq


def write_file(seq):
    """Here a fasta output gets written

    :param seq: string; containing an dna sequence
    :return: none

    Here 50 nucleotides get printed on one line.
    """
    with open(sys.argv[4], 'w') as file:
        file.write(f'>{sys.argv[7]}\n')
        for i in range(0, len(seq), 50):
            file.write(f'{seq[i:i + 50]}\n')

def main():
    """this is the main function of the script"""
    if os.path.exists(f'{sys.argv[5]}'):
        os.remove(f'{sys.argv[5]}')
    log = open(f'{sys.argv[5]}', 'a')
    sys.stdout = log
    fasta = open(sys.argv[1]).read().upper().strip().split('>')

    depth_file = open(sys.argv[6]).read().strip().split("\n")

    dic = parse_depth_file(depth_file[1:])

    seq = list(parse_fasta(fasta))
    seq = replace_low_depth_nuc(dic, seq)

    positions = open(sys.argv[2]).read().strip().split('\n')

    if re.match('\d+_[A-Z]_[A-Z]', positions[0]):
        sequence = alter_sequence(seq, positions)
    else:
        sequence = alter_sequence(seq, positions[1:])
    start_stop = sys.argv[3].split('-')
    sequence = ''.join(sequence[(int(start_stop[0]) - 1):int(start_stop[1])])
    write_file(sequence)


if __name__ == "__main__":
    main()
