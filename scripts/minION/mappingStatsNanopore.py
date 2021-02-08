import argparse
import sys

import pysam


def readSamAndGenerateStats(in_sam, fo, outputDetail):
    dicMap = {}  # For each read length the number of mapped reads
    dicUMap = {}  # For each read length the number of unmapped reads
    dicSec = {}  # For each read length the number of reads with
    # at least one secondary
    dicChim = {}  # For each read length the number of reads where
    # the primary alignment is chimeric
    dicCov = {}  # For each read length the sum of coverage on the reference
    dicReadsIDSecOrChim = {}  # In case the bam is not readname sorted will
    # be used to correct the dicChim and dicSec
    nReads = 0  # Total number of reads
    needToReReadTheBAM = False  # This is used if the BAM isn't readname sorted
    with pysam.Samfile(in_sam, 'r') as f:
        # Initialisation of variables
        currentReadName = ""
        isMapped = False
        length = 0
        hasSecondaryAlignment = False
        # mapq = 0
        # diffMapq = 0
        primaryAisChimeric = False
        mappedLength = 0
        reportValues = False
        for read in f.fetch():
            readname = read.qname.split("/")[0]
            if readname != currentReadName:
                # This is a new alignment
                if reportValues:
                    # You store the result or the currentReadName
                    if currentReadName in dicReadsIDSecOrChim:
                        # If the BAM is not readname sorted
                        # Then dicSec and dicChim will be updated after
                        hasSecondaryAlignment = False
                        primaryAisChimeric = False
                    dicMap[length] = dicMap.setdefault(length, 0) \
                        + int(isMapped)
                    dicUMap[length] = dicUMap.setdefault(length, 0) \
                        + int(not(isMapped))
                    dicSec[length] = dicSec.setdefault(length, 0) \
                        + int(hasSecondaryAlignment)
                    dicChim[length] = dicChim.setdefault(length, 0) \
                        + int(primaryAisChimeric)
                    dicCov[length] = dicCov.setdefault(length, 0) \
                        + mappedLength
                # Update the values for the new read
                currentReadName = readname
                isMapped = not(read.is_unmapped)
                if not(read.is_secondary) and not(read.is_supplementary):
                    # This is what will contribute to the coverage
                    # This is what happens always when
                    # the BAM is readname sorted.
                    length = read.query_length
                    nReads += 1
                    reportValues = True
                    if isMapped:
                        mappedLength = read.reference_length
                        if readname in dicReadsIDSecOrChim:
                            dicReadsIDSecOrChim[readname]["length"] = length
                    else:
                        mappedLength = 0
                    hasSecondaryAlignment = False
                    # mapq = 0
                    # diffMapq = 0
                    primaryAisChimeric = False
                else:
                    # The BAM is not readname sorted,
                    # so the alignment is stored to
                    # adjust the dictionaries values at the end.
                    reportValues = False
                    if not(readname in dicReadsIDSecOrChim):
                        dicReadsIDSecOrChim[readname] = \
                            {"length": 0,
                             "hasChim": False,
                             "hasSec": False}
                    if read.is_secondary:
                        dicReadsIDSecOrChim[readname]["hasSec"] = True
                    elif read.is_supplementary:
                        dicReadsIDSecOrChim[readname]["hasChim"] = True
            else:
                if not(read.is_secondary) and not(read.is_supplementary):
                    # This may happen when the BAM is not readname sorted
                    # This is what will contribute to the coverage
                    length = read.query_length
                    nReads += 1
                    reportValues = True
                    if isMapped:
                        mappedLength = read.reference_length
                        if readname in dicReadsIDSecOrChim:
                            dicReadsIDSecOrChim[readname]["length"] = length
                    else:
                        mappedLength = 0
                    hasSecondaryAlignment = False
                    # mapq = 0
                    # diffMapq = 0
                    primaryAisChimeric = False
                elif read.is_secondary:
                    hasSecondaryAlignment = True
                elif read.is_supplementary:
                    # The primaryAisChimeric is updated only when the BAM is
                    # readname sorted
                    primaryAisChimeric = True
        if reportValues:
            # You store the result or the currentReadName
            if currentReadName in dicReadsIDSecOrChim:
                # If the BAM is not readname sorted
                # Then dicSec and dicChim will be updated after
                hasSecondaryAlignment = False
                primaryAisChimeric = False
            dicMap[length] = dicMap.setdefault(length, 0) \
                + int(isMapped)
            dicUMap[length] = dicUMap.setdefault(length, 0) \
                + int(not(isMapped))
            dicSec[length] = dicSec.setdefault(length, 0) \
                + int(hasSecondaryAlignment)
            dicChim[length] = dicChim.setdefault(length, 0) \
                + int(primaryAisChimeric)
            dicCov[length] = dicCov.setdefault(length, 0) \
                + mappedLength
    if len(dicReadsIDSecOrChim) > 0:
        # the BAM is not readname sorted
        sys.stderr.write("The SAM file was not qname sorted.\n")
        if 0 in [d["length"] for d in dicReadsIDSecOrChim.values()]:
            sys.stderr.write("Some secondaries or chimerics were found after")
            sys.stderr.write(" the primary.\nThe file will be read")
            sys.stderr.write(" for a second time.\n")
            needToReReadTheBAM = True
        elif (sum(dicSec)+sum(dicChim)) > 0:
            # It happened at least once that their was a secondary/chimeric
            # following the primary.
            # First we check the Chimeric
            lengthsWithNewChim = [d["length"]
                                  for d in dicReadsIDSecOrChim.values()
                                  if d["hasChim"]]
            lengthsWithChim = [l for l in dicChim if dicChim[l] > 0]
            overlap = [i for i in lengthsWithChim if i in lengthsWithNewChim]
            if len(overlap) > 0:
                numberOfMappedReadsPerLength = [dicMap[l] for l in overlap]
                if max(numberOfMappedReadsPerLength) > 1:
                    # You will not be able to determine if they are the same
                    # reads or from different reads
                    sys.stderr.write("Some chimerics origin was ambigous.")
                    sys.stderr.write("\nThe file will be read for a")
                    sys.stderr.write(" second time.\n")
                    needToReReadTheBAM = True
                else:
                    for rn in dicReadsIDSecOrChim:
                        if dicReadsIDSecOrChim[rn]["length"] in overlap:
                            # The chimeric reads hare already included
                            # in dicChim
                            # Then the dicReadsIDSecOrChim is adjusted
                            dicReadsIDSecOrChim[rn]["hasChim"] = False
            if not needToReReadTheBAM:
                # We now check the Secondaries:
                lengthsWithNewSec = [d["length"]
                                     for d in dicReadsIDSecOrChim.values()
                                     if d["hasSec"]]
                lengthsWithSec = [l for l in dicSec if dicSec[l] > 0]
                overlap = [i for i in lengthsWithSec
                           if i in lengthsWithNewSec]
                if len(overlap) > 0:
                    numberOfMappedReadsPerLength = [dicMap[l] for l in overlap]
                    if max(numberOfMappedReadsPerLength) > 1:
                        # You will not be able to determine
                        # if they are the same
                        # reads or from different reads
                        sys.stderr.write("Some secondaries origin was ")
                        sys.stderr.write("ambigous.\nThe file will be read")
                        sys.stderr.write(" for a second time.\n")
                        needToReReadTheBAM = True
                    else:
                        for rn in dicReadsIDSecOrChim:
                            if dicReadsIDSecOrChim[rn]["length"] in overlap:
                                # The secondary reads hare already included
                                # in dicSec
                                # Then the dicReadsIDSecOrChim is adjusted
                                dicReadsIDSecOrChim[rn]["hasSec"] = False
        if not needToReReadTheBAM:
            # The dictionnaries just need to be adjusted:
            for rn in dicReadsIDSecOrChim:
                length = dicReadsIDSecOrChim[rn]["length"]
                dicSec[length] = dicSec.setdefault(length, 0) \
                    + int(dicReadsIDSecOrChim[rn]["hasSec"])
                dicChim[length] = dicChim.setdefault(length, 0) \
                    + int(dicReadsIDSecOrChim[rn]["hasChim"])
    if needToReReadTheBAM:
        with pysam.Samfile(in_sam, 'r') as f:
            for read in f.fetch():
                readname = read.qname.split("/")[0]
                if not(readname in dicReadsIDSecOrChim):
                    # I need to put each read in the dic because
                    # I need to have the real read length which
                    # is only present in the primary
                    dicReadsIDSecOrChim[readname] = \
                        {"length": 0,
                         "hasChim": False,
                         "hasSec": False}
                if read.is_secondary:
                    dicReadsIDSecOrChim[readname]["hasSec"] = True
                elif read.is_supplementary:
                    dicReadsIDSecOrChim[readname]["hasChim"] = True
                else:
                    dicReadsIDSecOrChim[readname]["length"] = read.query_length
        dicSec = {}
        dicChim = {}
        for rn in dicReadsIDSecOrChim:
            length = dicReadsIDSecOrChim[rn]["length"]
            dicSec[length] = dicSec.setdefault(length, 0) \
                + int(dicReadsIDSecOrChim[rn]["hasSec"])
            dicChim[length] = dicChim.setdefault(length, 0) \
                + int(dicReadsIDSecOrChim[rn]["hasChim"])
    fo.write("%i\tTotal reads\n" % nReads)
    totMappedReads = sum(dicMap.values())
    fo.write("%i\tMapped reads (%.2f %%)\n"
             % (totMappedReads, float(totMappedReads)/nReads*100))
    totUniquelyMappedReads = totMappedReads - sum(dicSec.values())
    fo.write("%i\tUniquely mapped reads (%.2f %%) \
the secondary alignment may be MAPQ 0\n"
             % (totUniquelyMappedReads,
                float(totUniquelyMappedReads)/nReads*100))
    totPrimaryChimer = sum(dicChim.values())
    fo.write("%i\tMapped reads as chimers (%.2f %%)\n"
             % (totPrimaryChimer, float(totPrimaryChimer)/nReads*100))
    nBases = sum([k*dicMap[k] for k in dicMap]) \
        + sum([k*dicUMap[k] for k in dicUMap])
    nBasesCovered = sum(dicCov.values())
    fo.write("%i\tTotal bases sequenced\n" % nBases)
    fo.write("%i\treference bases with coverage (%.2f %%)\n"
             % (nBasesCovered, float(nBasesCovered)/nBases*100))
    fo.write("Are counted here only the coverage generated by primary  alignment\
and only the best quality alignment in case of chimer.\n\
In addition, this coverage is larger than the real coverage\
which remove indels from the read.\n")
    if outputDetail != "":
        with open(outputDetail, 'w') as fo2:
            fo2.write("length\tnMapped\tnUnmapped\t\
nSecondaryAlignment\tnPrimaryIsChimeric\tCoverage\n")
            for l in sorted(set([k for k in dicMap] + [k for k in dicUMap])):
                fo2.write("%i\t%i\t%i\t%i\t%i\t%i\n"
                          % (l, dicMap.setdefault(l, 0),
                             dicUMap.setdefault(l, 0),
                             dicSec.setdefault(l, 0),
                             dicChim.setdefault(l, 0),
                             dicCov.setdefault(l, 0)))


argp = argparse.ArgumentParser(
  description=("Generate mapping stats from minimap2 sam"))
argp.add_argument('sam', default=None,
                  help="Input sam from minimap2.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output mapping report.")
argp.add_argument('--outputDetail', default="",
                  help=("Detailed output per size of reads."))
args = argp.parse_args()
readSamAndGenerateStats(args.sam, args.output, args.outputDetail)
