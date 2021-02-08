import argparse
import sys
import gzip


def store_in_dic(input_file, dic={}, norm_factor=1):
    if input_file.endswith('gz'):
        f = gzip.open(input_file, 'rb')
    else:
        f = open(input_file, 'r')
    for line in f:
        if isinstance(line, bytes):
            line = str(line, 'utf-8')
        ls = line.split('\t')
        possible_bp_full_pos = ls[0]
        possible_bp_seq = ls[1]
        dic[possible_bp_full_pos + '::' + possible_bp_seq] = \
            dic.setdefault(possible_bp_full_pos + '::' + possible_bp_seq, 0) + int(ls[3])/norm_factor
    return(dic)

def generate_fasta(dic, fo):
    for id, cov in sorted(dic.items(), key=lambda x: x[1], reverse=True):
        sequence = id.split('::')[-1]
        fo.write('>{}::{}\n'.format(id, cov))
        fo.write('{}\n'.format(sequence))

argp = argparse.ArgumentParser(
    description=("Combine BP files to generate fasta sorted by decrease coverage"))
argp.add_argument('--input', default=None, nargs='+',
                  required=True,
                  help="Input report(s) can be gzip")
argp.add_argument('--normFactors', default=None, nargs='+',
                  help="Space separated factors to use to normalize the coverage for each input (default is 1)")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output report.")
args = argp.parse_args()

if args.normFactors is None:
    norm_factors = [1] * len(args.input)
else:
    assert len(args.normFactors) == len(args.input), \
        "The number of normFactors must be the same as the number of input files."
dic = {}
for i, report in enumerate(args.input):
    dic = store_in_dic(report, dic, norm_factors[i])

generate_fasta(dic, args.output)
