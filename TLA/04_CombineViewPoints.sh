#!/bin/bash -l

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name postTLA
#SBATCH --chdir /scratch/ldelisle/TLA

# This script will combine information from all viewpoints

gitHubDirectory=$1

path="$PWD/"
pathForFastq=${path}fastq/
pathForSRATable=${gitHubDirectory}/TLA/sraTable.txt
pathForScripts=${gitHubDirectory}/scripts/TLA/

genome=mm10

module purge
module load gcc/7.4.0 # required for python
module load python/3.7.3

samples=$(cat $pathForSRATable | awk '{print $2}')

# Put here comma separated backBones
backBones=TgN3840_fosmid
backBonesSpace=$(echo $backBones | tr "," " ")

currentFastqSuffix="_unaligned_CATG_split_no${genome}"
for backBone in $backBonesSpace; do
  currentFastqSuffix=${currentFastqSuffix}no${backBone}
done

bpfiles=""
for sample in $samples; do
  currentFastq=${sample}${currentFastqSuffix}
  for backBone in $backBonesSpace $genome; do
    what="${backBone}"
    if [ "$backBone" = "$genome" ]; then
      what=$genome
    fi
    fullFile=${currentFastq}_findBP_full_on${what}.txt
    bpfiles="$bpfiles ${sample}/04_rmMotif_${fullFile}"
  done
done

# We compile all breakpoint files
# to generate fasta sorted by decrease coverage
python ${pathForScripts}generate_fasta.py --input $bpfiles --output BP_fromAll.fa

# We will compute the coverage of the top 50 breakpoints
# in each fastq
head -n 100 BP_fromAll.fa > BP_fromAll_top50.fa
input=BP_fromAll_top50
for sample in $samples; do
  python ${pathForScripts}reportOnBP.py --files ${pathForFastq}/${sample}.fastq.gz \
    --bpfile ${input}.fa --output ${input}_${sample}.txt
  if [ ! -e ${input}_all.txt ]; then
      cat ${input}_${sample}.txt | awk -v OFS="\t" -v name=$sample 'BEGIN{print "bp", name"_fwd_uniq", name"_fwd_tot", name"_rev_uniq", name"_rev_tot"}{print $1, $2, $3, $4, $5}' > ${input}_all.txt
  else
      TMPFILE=$(mktemp)
      cat ${input}_${sample}.txt | awk -v OFS="\t" -v name=$sample 'BEGIN{print name"_fwd_uniq", name"_fwd_tot", name"_rev_uniq", name"_rev_tot"}{print $2, $3, $4, $5}' | paste ${input}_all.txt - > $TMPFILE
      mv $TMPFILE ${input}_all.txt
  fi
done
# TODO

# Out of this BP_fromAll_top50_all.txt
# We select BP which are found in both samples
# Which does not have a strand bias (more than 10 times)
# We end up with this list:
# TgN3840:44996__to__::TCGCTCACTCGAGGGCTTCGCCCT::2562.0 = __to__TgN3840:0::TCGCTCACTCGAGGGCTTCGCCCT::6865.0
# __to__TgN3840:0::TCCTTACCTCGAGGGCTTCGCCCT::5285.0
# chr10:97019824__to__::AAATCCTTACCTCGAGGGCTTCGC::5118.0
# chr2:75138826__to__::AAAGAATTTATTGAGTGTTGCTTT::668.0 = TgN3840:18174__to__::AAAGAATTTATTGAGTGTTGCTTT::668.0
# __to__chr10:97019221::ATTCAACAGTTGTTGAAATAATGT::436.0

# The next one are likely false positive:
# __to__TgN3840:39999::GCTTTTTGGCCTCTGTCGTTTCCT::231.0
# chr2:181919135__to__::TATGTATGTATGTATGTATGTATG::149.0
# TgN3840:20891__to__::TGTGTGTGTGTGAGAGAGAGAGAG::144.0
# TgN3840:28813__to__::AAAAAAAAAAAGCCCAGGACCAGA::128.0
# __to__TgN3840:20837::TGTGTGTGTGTGTGTGTGTGTGTG::125.0
# __to__TgN3840:28813::GAAAGAAAGAAAGAAAGAAAGAAA::125.0

# Except this one which is the complement of __to__chr10:97019221
# TgN3840:26606__to__::CATAATTAAATTCAACAGTTGTTG::115.0


# For the 6 BP we can check the consensus:
# We need to be careful that what is outside of CATG CATG is probably wrong.

# For the tail to head of the fosmid:
# grep "TgN3840:44996__to__" */consensus*
# TgN3840_CS38/consensusOn04_TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:TgN3840:44996__to__    TCGCTCACTCGAGGGCTTCGCCCT        43      103     GATAACCCCAAGGGAAGTTTTTTCAGGCATCGTGTGTAAGCAGAATATATAAGTGCTGTTCCCTGGTGCTTCCTCGCTCACTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGG
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:TgN3840:44996__to__        TCGCTCACTCGAGGGCTTCGCCCT        305     2459    GTGGATAACCCCAAGGGAAGTTTTTTCAGGCATCGTGTGTAAGCAGAATATATAAGTGCTGTTCCCTGGTGCTTCCTCGCTCACTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGGAC


# For the left breakpoint:
# grep "__to__TgN3840:0" */consensus* | grep TCCTTACCTCGAGGGCTTCGCCCT
# TgN3840_CS38/consensusOn04_TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:__to__TgN3840:0TCCTTACCTCGAGGGCTTCGCCCT        114     724     CCATGTAGAAAATTGCAAGTAAATCCTTACCTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGGACAGACCACATCATG
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:__to__TgN3840:0    TCCTTACCTCGAGGGCTTCGCCCT        358     4561    ACATACGTTCCGCCATTCCTATGCGATGCACATGTAGAAAATTGCAAGTAAATCCTTACCTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGGACAGACCACATCATGTTATAT
# grep "chr10:97019824__to__" */consensus*
# TgN3840_CS38/consensusOn04_TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:chr10:97019824__to__ AAATCCTTACCTCGAGGGCTTCGC    80      552     CCATGTAGAAAATTGCAAGTAAATCCTTACCTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGGACAGACCACATCATG
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:chr10:97019824__to__AAATCCTTACCTCGAGGGCTTCGC 282     4566    TCGGGACACCACATACGTTCCGCCATTCCTATGCGATGCACATGTAGAAAATTGCAAGTAAATCCTTACCTCGAGGGCTTCGCCCTGTCGCTCGACTGCGGCGAGCACTACTGGCTGTAAAAGGACAGACCACATCATG

# For the right breakpoint:
# grep "chr2:75138826__to__" */consensus*
# TgN3840_CS38/consensusOn04_TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:chr2:75138826__to__  AAAGAATTTATTGAGTGTTGCTTT    178     600     AGCATTACGGCCCATGTGACTTTCCGCGGTCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACT
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:chr2:75138826__to__ AAAGAATTTATTGAGTGTTGCTTT 44      68      CCATGTGACTTTCCGCGGTCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACT
# grep "__to__chr10:97019221" */consensus*
# TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:__to__chr10:97019221 ATTCAACAGTTGTTGAAATAATGT    67      331     GTCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACTCGCGTGGATTCAGTGTGGTTTTCTAGCATACATATGATCACACAATGA
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onmm10.txt:__to__chr10:97019221ATTCAACAGTTGTTGAAATAATGT 40      105     TCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACTCGCGTGGATTCAGTGTGGTTTTCTAGCATACATATGATCACACAATG
# grep "TgN3840:26606__to__" */consensus*
# TgN3840_CS38/consensusOn04_TgN3840_CS38_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:TgN3840:26606__to__    CATAATTAAATTCAACAGTTGTTG        41      87      TTGGATTTTCCCGGGTCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACTCGCGTGGATTCAGTGTGGTTTTCTAGCATACATA
# TgN3840_vector/consensusOn04_TgN3840_vector_unaligned_CATG_split_nomm10noTgN3840_fosmid_findBP_full_onTgN3840_fosmid.txt:TgN3840:26606__to__        CATAATTAAATTCAACAGTTGTTG        19      28      TCTCAATAAAAAAAAAAGAATTTATTGAGTGTTGCTTTCATCATAATTAAATTCAACAGTTGTTGAAATAATGTGTAAATGAACTCAGACTCGCGTGGATTCAGTGTGGTTTTCTAGCATACA
