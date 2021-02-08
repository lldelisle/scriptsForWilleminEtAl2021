# These analysis were done on our galaxy server.
# Here are the command lines and version of tools:

# ChIPmentation and ChIP:

# cutadapt version 1.16
cutadapt  -j ${GALAXY_SLOTS:-1} -a 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' --output='out1.fq' \
  --error-rate=0.1 --times=1 --overlap=3 --minimum-length=15   --quality-cutoff=30 'input.fq'  > report.txt

# bowtie2 2.3.4.1 (the reference is either mm10 or the mutant genome)
bowtie2  -p ${GALAXY_SLOTS:-4}  -x 'mm10' -U 'input_f.fastq' 2> 'stats.txt'  | samtools sort -@${GALAXY_SLOTS:-2} \
  -O bam -o 'out.bam'

# When mapping to mm10:
samtools view -o 'output.bam' -h   -b  -q 30 input.bam

# macs2 version 2.1.1.20160309 
macs2 callpeak  --name 'MACS2'   -t 'input.bam' --format BAM --gsize '1870000000' --keep-dup '1' --bdg \
  --qvalue '0.05' --mfold '5' '50' --bw '300'

# Merging replicates:
# First scale:
awk -v s=scale 'BEGIN {OFS="\t"}{if (substr($0,0,5)=="track"){print $0}else{$4=$4/s;print$0}}' input.bdg > output.bdg
# Then sort:
sortBed -i 'input.bdg'  > 'output.bdg'
# Then do the mean:
unionBedGraphs  -filler '0'  -i 'rep1.bdg' 'rep2.bdg' | awk '{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > 'output.bdg'
