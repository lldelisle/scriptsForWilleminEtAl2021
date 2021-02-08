mutantGenomeDirectory="/scratch/ldelisle/mutantGenome/"
pathForChIPGEO="/scratch/ldelisle/ChIP_GEO/"
gitHubDirectory="/home/ldelisle/softwares/scriptsForWilleminEtAl2021/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10
pathForWtFasta=${pathForFasta}${genome}.fa

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.6"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.6 hicexplorer=3.6 pygenometracks=3.6
fi
conda activate hicexplorer3.6

# The peaks are extended 100bp each direction
sample1=E12_Limbs_TgN3840_CTCF_rep1
cat ${pathForChIPGEO}/${sample1}_on_TgN3840.bed | awk -v OFS="\t" '$1=="chr10" && $2 < 97880000 && $3 > 95480001 {print $1, $2 - 100, $3 + 100, $4}' > ${sample1}_on_TgN3840_plotted_extended.bed
# The fasta is extracted
bedtools getfasta -fi ${mutantGenomeDirectory}/TgN3840/mm10_TgN3840.fa -bed ${sample1}_on_TgN3840_plotted_extended.bed > ${sample1}_on_TgN3840_plotted_extended.fa
# On the website of http://insulatordb.uthsc.edu/ in CTCFBS Prediction Tool the fasta is uploaded.
# The table output for MIT_LM7 is copied into a file:
OUTPUT=${gitHubDirectory}/CTCF/${sample1}_on_TgN3840_plotted_extended_insdb_output.txt
# The orientation of the CTCF motif is deduced and put to . if score is negative
cat $OUTPUT | awk -v OFS="\t" '{split($3, a, ":|-"); if($7 < 0){$6="."}; print a[1], a[2] + $4 + $5/2, a[2] + $4 + $5/2 + 1, $3, $7, $6}' > ${gitHubDirectory}/CTCF/${sample1}_on_TgN3840.bed
sample2=E12_Limbs_TgN3840_CTCF_rep2
cat ${pathForChIPGEO}/${sample2}_on_TgN3840.bed | awk -v OFS="\t" '$1=="chr10" && $2 < 97880000 && $3 > 95480001 {print $1, $2 - 100, $3 + 100, $4}' > ${sample2}_on_TgN3840_plotted_extended.bed
# The fasta is extracted
bedtools getfasta -fi ${mutantGenomeDirectory}/TgN3840/mm10_TgN3840.fa -bed ${sample2}_on_TgN3840_plotted_extended.bed > ${sample2}_on_TgN3840_plotted_extended.fa
# On the website of http://insulatordb.uthsc.edu/ in CTCFBS Prediction Tool the fasta is uploaded.
# The table output for MIT_LM7 is copied into a file:
OUTPUT=${gitHubDirectory}/CTCF/${sample2}_on_TgN3840_plotted_extended_insdb_output.txt
# The orientation of the CTCF motif is deduced and put to . if score is negative
cat $OUTPUT | awk -v OFS="\t" '{split($3, a, ":|-"); if($7 < 0){$6="."}; print a[1], a[2] + $4 + $5/2, a[2] + $4 + $5/2 + 1, $3, $7, $6}' > ${gitHubDirectory}/CTCF/${sample2}_on_TgN3840.bed

cat ${gitHubDirectory}/CTCF/${sample1}_on_TgN3840.bed ${gitHubDirectory}/CTCF/${sample2}_on_TgN3840.bed | sort -k1,1  -k2,2n | awk '
{
    if(cur_chr == $1 && cur_start == $2 && cur_end == $3 && cur_strand == $6){
        if(cur_score < $5){
            bestLine = $0
        }
    } else {
        if (NR > 1){
            print bestLine
        }
        bestLine=$0
        cur_chr=$1
        cur_start=$2
        cur_end=$3
        cur_strand=$6
    }
}' > ${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840.bed

# Create a bed with rgb field corresponding to motif orientation of CTCF:
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" -v colorOther="0,0,0" '
{
  if ($6 == "+"){
    color = colorPos
  } else if ($6 == "-" ){
    color = colorNeg
  } else {
    color = colorOther
  }
  print $1, $2, $3, $4, $5, $6, $2, $2, color
}' ${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840.bed > ${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840_colored.bed

# For the wt:
# The peaks are extended 100bp each direction
sample=E12_Limbs_Wt_CTCF
cat ${pathForChIPGEO}/${sample}.bed | awk -v OFS="\t" '($1=="chr10" && $2 < 97880000 && $3 > 95480001)||($1=="chr2" && $2 < 75800000 && $3 > 73640001) {print $1, $2 - 100, $3 + 100, $4}' > ${sample}_plotted_extended.bed
# The fasta is extracted
bedtools getfasta -fi ${pathForWtFasta} -bed ${sample}_plotted_extended.bed > ${sample}_plotted_extended.fa
# On the website of http://insulatordb.uthsc.edu/ in CTCFBS Prediction Tool the fasta is uploaded.
# The table output for MIT_LM7 is copied into a file:
OUTPUT=${gitHubDirectory}/CTCF/${sample}_plotted_extended_insdb_output.txt
# The orientation of the CTCF motif is deduced and put to . if score is negative
cat $OUTPUT | awk -v OFS="\t" '{split($3, a, ":|-"); if($7 < 0){$6="."}; print a[1], a[2] + $4 + $5/2, a[2] + $4 + $5/2 + 1, $3, $7, $6}' > ${gitHubDirectory}/CTCF/${sample}.bed

# Create a bed with rgb field corresponding to motif orientation of CTCF:
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" -v colorOther="0,0,0" '
{
  if ($6 == "+"){
    color = colorPos
  } else if ($6 == "-" ){
    color = colorNeg
  } else {
    color = colorOther
  }
  print $1, $2, $3, $4, $5, $6, $2, $2, color
}' ${gitHubDirectory}/CTCF/${sample}.bed > ${gitHubDirectory}/CTCF/${sample}_colored.bed

