import argparse
import sys


def checkBP(input_file, fo, bpnames, bpseqs):
    with open(input_file, 'r') as f:
        for line in f:
            ls = line.strip().split('\t')
            seqBP = ls[1]
            if bpnames is not None and ls[0] not in bpnames:
                continue
            if bpseqs is not None and seqBP not in bpseqs:
                continue
            seqDic = eval(ls[-1])
            seqDicOrder = {}
            for seq in seqDic:
                iBP = seq.find(seqBP)
                seqDicOrder.setdefault(iBP, []).append((seq, seqDic[seq]))
            allIBP = [iBP for iBP in seqDicOrder]
            maxI = max(allIBP)
            fo.write('\t'.join(ls[:-1]) + '\n')
            for iBP in sorted(allIBP, reverse=True):
                for seq, n in seqDicOrder[iBP]:
                    fo.write('{}\t{}\n'.format(n, ' ' * (maxI - iBP) + seq))
            fo.write('\n')

argp = argparse.ArgumentParser(
  description=("Display all unique reads which lead to BP identification correctly aligned from full report"))
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
args = argp.parse_args()

checkBP(args.input, args.output, args.bpname, args.bpseq)
