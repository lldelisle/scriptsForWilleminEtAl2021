import argparse
import sys
import gzip

def one_basepair_appart(seq1, seq2):
    # First test is aligned but 1 mm:
    if len([base for i, base in enumerate(seq1) if base != seq2[i]]) == 1:
        return(True)
    # Second test is seq1 shifted 1bp to left:
    if len([base for i, base in enumerate(seq1[1:]) if base != seq2[i]]) == 0:
        return(True)
    # Third test is seq1 shifted 1bp to right:
    if len([base for i, base in enumerate(seq1[:-1]) if base != seq2[i + 1]]) == 0:
        return(True)
    
    return(False)

def filter_findBP_file_with_redundancy(inputFile, fo):
    # We assume the file is already sorted and
    # Secondary hits will be discarded
    written_bp = {}
    if inputFile.endswith('gz'):
        f = gzip.open(inputFile, 'rb')
    else:
        f = open(inputFile, 'r')
    for line in f:
        to_write = True
        if isinstance(line, bytes):
            line = str(line, 'utf-8')
        ls = line.split('\t')
        possible_bp_full_pos = ls[0]
        match_right = possible_bp_full_pos.startswith('__to__')
        chrom, pos = possible_bp_full_pos.strip('__to__').split(':')
        pos = int(pos)
        possible_bp_seq = ls[1]
        possible_matching_bp = ['{}:{}:{}'.format(chrom, p, match_right) for p in range(pos - 1, pos + 2)]
        matching_bp = [mbp for mbp in possible_matching_bp if mbp in written_bp]
        if len(matching_bp) > 0:
            for mbp in matching_bp:
                if one_basepair_appart(possible_bp_seq, written_bp[mbp]):
                    to_write = False
        if to_write:
            fo.write(line)
            written_bp['{}:{}:{}'.format(chrom, pos, match_right)] = possible_bp_seq
        
        

argp = argparse.ArgumentParser(
  description=("Remove all BP which are already present upper in the file and has only one mismatch with it (or 1bp appart)."))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Sorted input file to filter (output of findBP sorted by sort -k4,4nr for example), can be gzip.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output filtered file.")
args = argp.parse_args()

filter_findBP_file_with_redundancy(args.input, args.output)
