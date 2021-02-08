# scriptsForWilleminEtAl2021

All the scripts needed to reproduce the NGS figures in Willemin et al. 2021 from raw data.

mutantGenome directory gives the scripts to build the mm10_TgN3840 mutant genome fasta file, bowtie2 indices and other input for the 4C-seq analysis.

4C directory contains a master bash script (pipeline.sh) with all the command lines to run the 4C-seq analysis (from fastq to bedgraph and bigwig).

HiC directory contains the slurm bash scripts which have been used to process the data from fastq to cool and TAD calling.

TLA directory contains the slurm bash scripts which have been used to process the data from fastq to coverage and breakpoint identification. It also contains bed files and fasta files of the fosmid.

minION_nCATS directory contains the slurm bash scripts which have been used to process the data from fast5 to statistics in supplementary table 4 and fasta files of selected reads. It also contains the output which have been used to generate the supplementary table 4 and the fasta file of the construct TLA-derived.

The ChIP and ChIPmentation analyses were done on our galaxyserver.

The CTCF directory contains a bash script which has been used to generate the bed file with motif orientation. It also contains the intermediate results: the bed with peaks in the region and the output of insulatorDB.

The plots directory contains a slurm bash script with all pyGenomeTracks commands and all other commands needed to compute statistics. It also contains a directory outputs with all figures before modification with illustrator.

Finally the scripts directory contains all scripts which were used in the analysis.
