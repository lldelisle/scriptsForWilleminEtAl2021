import argparse
import sys
import pysam

def findBP(inputBAM, minCov, minCovTot, fo, outputFull):
    basesAroundBP = 12
    possibleBP = {}
    bamfile = pysam.AlignmentFile(inputBAM, "rb")
    refLengths = bamfile.lengths
    for read in bamfile.fetch():
        first_cigar = read.cigartuples[0]
        # We only catch reads with at least 12S
        if first_cigar[0] == 4 and first_cigar[1] >= basesAroundBP:
            possible_bp_pos = read.reference_start
            possible_bp_ref_name = bamfile.get_reference_name(read.reference_id)
            if 'twice' in possible_bp_ref_name and possible_bp_pos > refLengths[read.reference_id] / 2:
              possible_bp_pos -= int(refLengths[read.reference_id] / 2)
            possible_bp_full_pos = '__to__' + possible_bp_ref_name + ':{}'.format(possible_bp_pos)
            bp_index = first_cigar[1]
            possible_bp_seq = read.query_sequence[bp_index - basesAroundBP:bp_index + basesAroundBP]
            possibleBP.setdefault(possible_bp_full_pos, {}).setdefault(possible_bp_seq, {})[read.query_sequence] = \
                possibleBP.setdefault(possible_bp_full_pos, {}).setdefault(possible_bp_seq, {}).setdefault(read.query_sequence, 0) + 1
        last_cigar = read.cigartuples[-1]
        # We only catch reads with at least 12S
        if last_cigar[0] == 4 and last_cigar[1] >= basesAroundBP:
            possible_bp_pos = read.reference_end
            possible_bp_ref_name = bamfile.get_reference_name(read.reference_id)
            if 'twice' in possible_bp_ref_name and possible_bp_pos > refLengths[read.reference_id] / 2:
              possible_bp_pos -= int(refLengths[read.reference_id] / 2)
            possible_bp_full_pos = possible_bp_ref_name + ':{}__to__'.format(possible_bp_pos)
            bp_index = read.query_length - last_cigar[1]
            possible_bp_seq = read.query_sequence[bp_index - basesAroundBP:bp_index + basesAroundBP]
            possibleBP.setdefault(possible_bp_full_pos, {}).setdefault(possible_bp_seq, {})[read.query_sequence] = \
                possibleBP.setdefault(possible_bp_full_pos, {}).setdefault(possible_bp_seq, {}).setdefault(read.query_sequence, 0) + 1
    if outputFull is not None:
        fof = open(outputFull, 'w')    
    for possible_bp_full_pos in possibleBP:
        for possible_bp_seq in possibleBP[possible_bp_full_pos]:
            indexDic = possibleBP[possible_bp_full_pos][possible_bp_seq]
            nb_indices = len(indexDic)
            tot_coverage = sum([indexDic[k] for k in indexDic])
            if nb_indices >= minCov and tot_coverage >= minCovTot:
                fo.write('{}\t{}\t{}\t{}\n'.format(possible_bp_full_pos, possible_bp_seq,
                                                   nb_indices, tot_coverage))
                if outputFull is not None:
                    fof.write('{}\t{}\t{}\t{}\t{}\n'
                              ''.format(possible_bp_full_pos, possible_bp_seq,
                                        nb_indices, tot_coverage,
                                        indexDic))


argp = argparse.ArgumentParser(
  description=("Find BP from bam file"))
argp.add_argument('--input', default=None,
                  required=True,
                  help="Input bam with index")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output report.")
argp.add_argument('--outputFull', default=None,
                  help="Output full report.")
argp.add_argument('--minCov', default=0,
                  type=int,
                  help="Minimum coverage (identical position on reads are considered as one) to report a break-point.")
argp.add_argument('--minCovTot', default=0,
                  type=int,
                  help="Minimum total coverage to report a break-point.")
args = argp.parse_args()

findBP(args.input, args.minCov, args.minCovTot, args.output, args.outputFull)
