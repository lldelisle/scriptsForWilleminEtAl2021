import argparse
import sys
import gzip

def all_motif_to_dic(all_motif_file):
    # We assume each chromosome is only once
    dic = {}
    if all_motif_file.endswith('gz'):
        f = gzip.open(all_motif_file, 'r')
    else:
        f = open(all_motif_file, 'r')
    for line in f:
        if isinstance(line, bytes):
            line = str(line, 'utf-8')
        ls = line.split('\t')
        chrom = ls[0]
        dic[chrom] = [int(p) for p in ls[1].split(', ')]
    f.close()
    return(dic)

def filter_findBP_file_with_motif(inputFile, motifDic, fo):
    motif = 'CATG'
    if inputFile.endswith('gz'):
        f = gzip.open(inputFile, 'rb')
    else:
        f = open(inputFile, 'r')
    for line in f:
        if isinstance(line, bytes):
            line = str(line, 'utf-8')
        ls = line.split('\t')
        possible_bp_full_pos = ls[0]
        chrom, pos = possible_bp_full_pos.strip('__to__').split(':')
        pos = int(pos)
        possible_bp_seq = ls[1]
        motif_around_bp = [i for i in range(pos - len(motif) - 1, pos + 1) if i in motifDic[chrom]]
        if len(motif_around_bp) > 0:
            # There is a motif before:
            # We need to count the number of matches with the motif
            # for example: __to__I7revfwd_twice:15546      CAGAGGTGACACGTGTAGCTGGAG
            # and motif_around_bp = [15543] so                      CATG
            # for example: I7revfwd_twice:16510__to__      AGAGGTGACGTGGGGCGTGTTAAT
            # and motif_around_bp = [16506] so                     CATG
            has_partial_match = False
            for motif_pos in motif_around_bp:
                starting_index = int(len(possible_bp_seq) / 2) + motif_pos - pos
                potential_motif = possible_bp_seq[starting_index: starting_index + len(motif)]
                matches = len([base for i, base in enumerate(potential_motif) if base == motif[i]])
                if matches == len(motif) - 1:
                    has_partial_match = True
            if not has_partial_match:
                fo.write(line)
        else:
            fo.write(line)
        

argp = argparse.ArgumentParser(
  description=("Remove all BP which are close to a motif in the reference and has only one mismatch with it."))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Input file to filter (output of findBP), can be gzip.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output filtered file.")
argp.add_argument('--motifFile', default=None,
                  help="File with all positions of motif one line per chromosome, output of findAllREPos_rep.py, can be gzip.")
args = argp.parse_args()

motifDic = all_motif_to_dic(args.motifFile)
filter_findBP_file_with_motif(args.input, motifDic, args.output)
