import sys
import argparse
from Bio import SeqIO
import re

def findAllMotifs(fasta_path, fo, motif_seq, verbose):
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        if verbose:
            print("Processing {}".format(seq_record.id))
        l = [m.start() for m in re.finditer(motif_seq, str(seq_record.seq).upper())]
        fo.write(seq_record.id + '\t')
        fo.write(str(l)[1:-1] + '\n')

argp = argparse.ArgumentParser(
  description=("Find all positions of a motif in a fasta file."))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Input fasta, can be gzip.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output positions (0-based, one line per chromosome).")
argp.add_argument('--motif', default=None, required=True,
                  help="Motif to look in the fasta")
argp.add_argument('--verbose', action='store_true',
                  help="Display messages to show the process.")
args = argp.parse_args()

findAllMotifs(args.input, args.output, args.motif, args.verbose)
