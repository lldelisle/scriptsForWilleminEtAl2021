import argparse
import sys
import gzip

def bedToDic(bedFile, dic={}):
    if bedFile.endswith('gz'):
        f = gzip.open(bedFile, 'r')
    else:
        f = open(bedFile, 'r')
    for line in f:
        if isinstance(line, bytes):
            line = str(line, 'utf-8')
        ls = line.split('\t')
        chrom = ls[0]
        if not chrom in dic:
            dic[chrom] = {}
        for i in range(int(ls[1]), int(ls[2])):
            dic[chrom][i] = True
    f.close()
    return(dic)

def filterFindBPFile(inputFile, bedDic, fo, polyAfilter):
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
        possible_bp_seq = ls[1]
        is_poly_A = (possible_bp_seq == 'A' * len(possible_bp_seq)) or \
            (possible_bp_seq == 'T' * len(possible_bp_seq))
        if int(pos) not in bedDic.setdefault(chrom, {}) and not (is_poly_A and polyAfilter):
            fo.write(line)
        


argp = argparse.ArgumentParser(
  description=("Remove all BP which overlap with a bed (for example repeat masker)."))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Input file to filter (output of findBP), can be gzip.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output filtered file.")
argp.add_argument('--polyA', action='store_true', help="Filter bp whose sequence is only A/T.")
argp.add_argument('--bed', default=None, nargs='+',
                  help="Bed(s) with intervals to remove (for example repeat masker), can be gzip.")
args = argp.parse_args()

bedDic = {}
if args.bed is not None:
    for bedFile in args.bed:
        bedDic = bedToDic(bedFile, bedDic)
filterFindBPFile(args.input, bedDic, args.output, args.polyA)
