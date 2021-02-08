import argparse
import sys
import gzip

def revComplement(s):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(s.upper()))

def fastaToDic(fastaFile):
    dic={}
    name=""
    with open(fastaFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
            elif len(line) == 0:
                next
            else:
                if name == "":
                    raise Exception("Invalid fasta file, sequence before name.")
                dic[line] = name
                if revComplement(line) != line:
                    dic[revComplement(line)] = 'rev_' + name
    return(dic)

def generateStats(filesList, fastaFile, fo, outputFull):
    dicWithSeq = fastaToDic(fastaFile)
    seqWithBP = {seq: {} for seq in dicWithSeq}
    for inputFile in filesList:
        with gzip.open(inputFile, 'rb') as f:
            iLine = 0
            for line in f:
                iLine += 1
                if iLine % 4 == 2:
                    lineSeq = str(line, 'utf-8').strip()
                    for seqBP in seqWithBP:
                        iSeqBP = lineSeq.find(seqBP)
                        if iSeqBP >= 0:
                            seqWithBP[seqBP].setdefault(iSeqBP, {})[lineSeq] = \
                                seqWithBP[seqBP].setdefault(iSeqBP, {}).setdefault(lineSeq, 0) + 1
    if outputFull is not None:
        fof = open(outputFull, 'w')
    bpNames = {dicWithSeq[k]: k for k in dicWithSeq}
    fwdBPNames = [bpName for bpName in bpNames if not bpName.startswith('rev_')]
    for bpName in fwdBPNames:
        seqBP = bpNames[bpName]
        difSeq = 0
        totSeq = 0
        for iSeqBP in sorted(seqWithBP[seqBP].keys()):
            for lineSeq in seqWithBP[seqBP][iSeqBP]:
                difSeq += 1
                totSeq += seqWithBP[seqBP][iSeqBP][lineSeq]
                if outputFull is not None:
                    fof.write(bpName + '\t{}\t{}\t{}\n'.format(iSeqBP,
                                                               seqWithBP[seqBP][iSeqBP][lineSeq],
                                                               lineSeq))
        fo.write(bpName + '\t{}\t{}'.format(difSeq, totSeq))
        if revComplement(seqBP) != seqBP:
            seqBP = bpNames['rev_' + bpName]
        difSeq = 0
        totSeq = 0
        for iSeqBP in sorted(seqWithBP[seqBP].keys()):
            for lineSeq in seqWithBP[seqBP][iSeqBP]:
                difSeq += 1
                totSeq += seqWithBP[seqBP][iSeqBP][lineSeq]
                if outputFull is not None:
                    fof.write('rev_' + bpName + '\t{}\t{}\t{}\n'.format(iSeqBP,
                                                                        seqWithBP[seqBP][iSeqBP][lineSeq],
                                                                        lineSeq))
        fo.write('\t{}\t{}\n'.format(difSeq, totSeq))


argp = argparse.ArgumentParser(
  description=("Generate stats from fastq.gz and fasta file"))
argp.add_argument('--files', default=None, nargs='+',
                  required=True,
                  help="Input fastq.gz.")
argp.add_argument('--bpfile', default=None, required=True,
                  help="Input fasta with break-points")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output report.")
argp.add_argument('--outputFull', default=None,
                  help="Output full report.")
args = argp.parse_args()

generateStats(args.files, args.bpfile, args.output, args.outputFull)
