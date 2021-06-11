import os
import pickle
import sys
import argparse

def getSummary(dic):
    totNumberOfReads = dic['total'] + dic['unmapped']
    mappedReads = dic['total']
    v = [ str(totNumberOfReads),
        str(mappedReads),
        "%.2f" % (float(mappedReads) / totNumberOfReads * 100) ]
    return(v)

def writeBigSummary(inputDir, fo):
    header = 'Name\tRead total\tRead mapped\tProp Read mapped\n'
    fo.write(header)
    for bamstatfile in os.listdir(inputDir):
        if '_filter_bamstat' in bamstatfile:
            name = bamstatfile[:-len('_filter_bamstat')]
            dic = pickle.load(open(os.path.join(inputDir, bamstatfile), 'rb'))
            summary = getSummary(dic)
            fo.write('\t'.join([name] + summary) + "\n")

argp=argparse.ArgumentParser(description='Write a table with mapping rate from a folder of results with bamstat files.')
argp.add_argument('--directory',default="./",help="A directory with files bamstat",type=str)
argp.add_argument('--output',default=sys.stdout,type=argparse.FileType('w'))

args = argp.parse_args()

writeBigSummary(args.directory, args.output)