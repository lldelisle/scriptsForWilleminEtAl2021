import sys
import os

import pyBigWig

path_4C_geo = sys.argv[1]
output_file_stats_4C = sys.argv[2]
plotted_region = (95480001, 97880000)

# Stats on 4C: nb of fragments covered, signal in trans etc.
with open(output_file_stats_4C, 'w') as fo:
    fo.write("file_name\tNb_frags\tTrans\tSignalInPlottedRegion\tFragmentsInPlottedRegion\n")
    # The segToFrag contains the scores stored as bw file normalized to million mapped reads
    for segToFrag in [f for f in os.listdir(path_4C_geo) if f.startswith("segToFrag") and f.endswith(".bw")]:
        bw = pyBigWig.open(os.path.join(path_4C_geo, segToFrag))
        all_count = 0
        all_signal = 0
        signal_chr10 = 0
        for chrom in bw.chroms():
            intervals = bw.intervals(chrom, 0, bw.chroms(chrom))
            all_count += len(intervals)
            signal = sum([v[2] for v in intervals])
            all_signal += signal
            if chrom == "chr10":
                signal_chr10 = signal
                signal_plotted_region = sum([v[2] for v in bw.intervals(chrom, *plotted_region)])
                frags_plotted_region = len(bw.intervals(chrom, *plotted_region))
        fo.write(f"{segToFrag}\t{all_count}\t{100 * (1 - signal_chr10 / all_signal)}\t{100 * signal_plotted_region / all_signal}\t{frags_plotted_region}\n")
