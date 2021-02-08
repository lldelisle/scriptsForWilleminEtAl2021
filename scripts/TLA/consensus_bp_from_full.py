import argparse
import sys


def consensusBP(input_file, fo, bpnames, bpseqs, removeDup, min_cov):
    with open(input_file, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            seqBP = ls[1]
            if bpnames is not None and ls[0] not in bpnames:
                continue
            if bpseqs is not None and seqBP not in bpseqs:
                continue
            seqDic = eval(ls[-1])
            consensusDic = {}
            for seq in seqDic:
                iBP = seq.find(seqBP)
                if removeDup:
                    n = 1
                else:
                    n = seqDic[seq]
                for i, l in enumerate(seq):
                    consensusDic.setdefault((i - iBP), {})[l] = \
                        consensusDic.setdefault((i - iBP), {}).setdefault(l, 0) + n
            allI = [i for i in consensusDic]
            fo.write('\t'.join(ls[:-1]) + '\t')
            for i in sorted(allI):
                most_frequent_letter = max(consensusDic[i], key=consensusDic[i].get)
                if consensusDic[i][most_frequent_letter] >= min_cov:
                    fo.write(most_frequent_letter)
            fo.write('\n')

argp = argparse.ArgumentParser(
  description=("Check BP from full report"))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Input full report")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output report.")
argp.add_argument('--bpname', default=None, nargs='+',
                  help="Name(s) of bp to check")
argp.add_argument('--bpseq', default=None, nargs='+',
                  help="Sequence(s) of bp to check")
argp.add_argument('--removeDup', action='store_true',
                  help="Consider identical reads as one.")
argp.add_argument('--minCoverage', default=0, type=int,
                  help="Minimum coverage to report in the consensus sequence.")
args = argp.parse_args()

consensusBP(args.input, args.output, args.bpname, args.bpseq, args.removeDup, args.minCoverage)
