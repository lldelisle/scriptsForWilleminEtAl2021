#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 8:00:00
#SBATCH --job-name plots
#SBATCH --chdir /scratch/ldelisle/plots/

gitHubDirectory=$1

path=$PWD/

pathForScripts=${gitHubDirectory}/scripts/
pathForTLA=/scratch/ldelisle/TLA/
pathForSRATableTLA=${gitHubDirectory}/TLA/sraTable.txt
samplesTLA=$(cat $pathForSRATableTLA | awk '{print $2}')

pathForHiC=/scratch/ldelisle/HiC/

pathFor4CGEO=/scratch/ldelisle/4C/analysis_TgN3840/toGEO/

pathForChIPGEO="/scratch/ldelisle/ChIP_GEO/"
pathForMutantGenomes="/scratch/ldelisle/mutantGenome/"

pathForFREECResults="/scratch/ldelisle/CNV/E12_Limbs_TgN3840_inputvscontrol_merged_neq8_chr2_71_78/"

pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10
pathForChromSizes=${pathForFasta}${genome}.fa.fai

pathMinION="/scratch/ldelisle/nCATS/TgN3840_nCATS/"

module purge
module load gcc/7.4.0
module load openblas/0.3.6-openmp 
module load r/3.6.0
module load samtools/1.9

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
exists=$(conda info --envs | awk '$1=="hicexplorer3.6"{print}' | wc -l)
if [ $exists -ne 1 ]; then
  conda create -y -n hicexplorer3.6 hicexplorer=3.6 pygenometracks=3.6
fi
conda activate hicexplorer3.6

# Integration borders:
echo -e "chr10\t97019221\t97019222\trightBP\nchr10\t97019824\t97019825\tleftBP" > integration_borders.bed
Rscript ${pathForScripts}/shiftAnnot_TgN3840.R integration_borders.bed 1 2 3 integration_borders_TgN3840_temp.bed
cat integration_borders_TgN3840_temp.bed | awk '$4=="leftBP"{if($2==$3){print $1"\t"$2"\t"$3+1"\t"$4}else{print}}' > integration_borders_TgN3840.bed

# Region between 2 new subTADs:
echo -e "chr10\t96120000\t97000000\tchr10\t97080000\t97400000\tinterSubTAD" > interSubTAD.bedpe

# Region used for fig2Cbottom:
echo -e "chr10\t97019018\t97019019\tfig2CbottomLeft\nchr10\t97020045\t97020046\tfig2CbottomRight" > fig2Cbottom.bed

# Find the restriction sites in fig2Cbottom
samtools faidx ~/genomes/fasta/mm10.fa chr10:97019018-97020046 | \
  sed 1d | tr -d "\n" | awk -v start=97019018 '
{
  IGNORECASE=1
  current_string=$1
  current_index=start - 1
  i=match(current_string,"CATG")
  while(i > 0){
    current_index = current_index + i - 1 + 4
    if ((substr(current_string, i -1, 1) == "A" ||
      substr(current_string, i -1, 1) == "G") &&
      (substr(current_string, i + 4, 1) == "C" ||
      substr(current_string, i + 4, 1) == "T")){
        print "chr10\t"current_index - 4 "\t" current_index + 1 "\tNlaIII+NspI"
    } else {
      print "chr10\t"current_index - 4 "\t" current_index"\tNlaII"
    }
    current_string=substr(current_string, i + 4)
    i=match(current_string,"CATG")
  }
}' > re_fig2Cbottom.bed

# Process genes:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
zcat gencode.vM23.annotation.gtf.gz | grep "protein_coding" > gencode.vM23.protcod.gtf
Rscript ${pathForScripts}/shiftAnnot_TgN3840.R gencode.vM23.protcod.gtf 1 4 5 gencode.vM23.protcod.mm10_TgN3840.gtf

# Shift the CTCF:
cat ${gitHubDirectory}/CTCF/E12_Limbs_Wt_CTCF_colored.bed | awk '$1=="chr10"{print}' > E12_Limbs_Wt_CTCF_onlychr10_colored.bed
Rscript ${pathForScripts}/shiftAnnot_TgN3840.R E12_Limbs_Wt_CTCF_onlychr10_colored.bed 1 2 3 6 E12_Limbs_Wt_CTCF_mm10_TgN3840_colored.bed

# Convert bedgraph to bigwig:
for f in ${pathForChIPGEO}/*bedgraph; do
  bigwig=${f/bedgraph/bw}
  if [[ $f = *"_on_TgN3840"* ]]; then
    my_sizes=${pathForMutantGenomes}/TgN3840/mm10_TgN3840.fa.fai
  else
    my_sizes=$pathForChromSizes
  fi
  if [ ! -e $bigwig ]; then
    tempfile=$(mktemp)
    cat $f | grep -v track | LC_COLLATE=C sort -k1,1 -k2,2n > $tempfile
    bedGraphToBigWig $tempfile $my_sizes $bigwig
  fi
done

# Concatenate peaks files:
for prot in CTCF RAD21; do
  sample=E12_Limbs_TgN3840_${prot}
  cat ${pathForChIPGEO}/${sample}_rep[12]_on_mm10.bed > ${pathForChIPGEO}/${sample}_rep1or2_on_mm10.bed
done

# Shift the H3K27Ac to mutant genome:
zcat ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_on_mm10.bedgraph.gz | awk '$1 == "chr10"{print}' > ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10.bedgraph
Rscript ${pathForScripts}/shiftAnnot_TgN3840.R ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10.bedgraph 1 2 3 ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10_shifted_Tg3840.bedgraph
my_sizes=${pathForMutantGenomes}/TgN3840/mm10_TgN3840.fa.fai
bedGraphToBigWig ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10_shifted_Tg3840.bedgraph $my_sizes ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10_shifted_Tg3840.bw

# Get bed for fosmid end in _TH_TT_HH
echo "track type=bed name=fosmid_boundaries_on_TH_TT_HH" > fosmid_boundaries_on_TH_TT_HH.bed
for i in {1..3}; do
  echo -e "TgN3840_fosmid_TH_TT_HH\t$((i * 44996))\t$((i * 44996 +1))" >> fosmid_boundaries_on_TH_TT_HH.bed
done

bedtools slop -i fosmid_boundaries_on_TH_TT_HH.bed -g ${gitHubDirectory}/TLA/fasta/TgN3840_fosmid_TH_TT_HH.fa.fai -b 500 > fosmid_boundaries_on_TH_TT_HH_1kb.bed

# For each figure the ini file is generated followed by the plot command
ini_file="fig1A-C.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 100000
scalebar_end_position = 75750000
height = 0.3
" > ${ini_file}

for genotype in Wt TgN3840; do
  sample=E12_Limbs_${genotype}_map_mm10
  echo "[Hi-C ${sample}]
file = ${pathForHiC}/${sample}/${sample}.40kb.cool
title = $sample 40kb
colormap = ['white', (1, 0.88, 0.66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
depth = 1040000
transform = no
show_masked_bins = true
rasterize = false
min_value = 0
max_value = 0.04
" >> ${ini_file}
  for size in 240 480; do
    if [ $size = "240" ]; then
      color=grey
    else
      color=black
    fi
    echo "[spacer]
height = 0.25

[${sample}.${size}kb_domains]
file = ${pathForHiC}/TADs_diff/${sample}.40kb.${size}kb_domains.bed
title = ${sample} 40kb ${size}kb_domains
color = ${color}
border_color = none
height = 0.5
labels = false
line_width = 1
display = interleaved
" >> ${ini_file}
  done
  echo "[spacer]
height = 1
" >> ${ini_file}
done
echo "[E12_Limbs_CTCF]
file = ${pathForChIPGEO}/E12_Limbs_Wt_CTCF_on_mm10.bw
title = Limbs_E12_Wt_CTCF
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = 6

[spacer]
height = 0.15

[E12_Limbs_CTCF_orientation]
file = ${gitHubDirectory}/CTCF/E12_Limbs_Wt_CTCF_colored.bed
title = MACS2_E12_Limbs_Wt_CTCF
height = 0.5
labels = false
line_width = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb

[spacer]
height = 1

[HoxD_TADs_minimal_transcripts]
file = ${pathForScripts}/annotations/HoxD_TADs_minimal_transcripts.gtf
title = HoxD TADs minimal transcripts
color = black
height = 0.5
labels = false
line_width = 0.5
display = collapsed

[spacer]
height = 1

[annotations_consensus_mm10_enhancers]
file = ${pathForScripts}/annotations/annotations_consensus_mm10_enhancers.bed
title = annotations_consensus_mm10_enhancers
color = black
height = 0.5
labels = true
line_width = 0.5
display = collapsed
file_type = bed

[spacer]
height = 0.5

[TgN3840_fosmid]
file = ${pathForScripts}/annotations/TgN3840_fosmid.bed
title = TgN3840_fosmid
color = black
height = 0.2
labels = true
line_width = 1
display = collapsed
file_type = bed
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:73640001-75800000 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 34.76 \
  --trackLabelFraction 0.3

ini_file="fig1D.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 5000
center = 75133816
height = 0.3

[E12_Limbs_CTCF]
file = ${pathForChIPGEO}/E12_Limbs_Wt_CTCF_on_mm10.bw
title = Limbs_E12_Wt_CTCF
height = 3
color = black
nans_to_zeros = true
number_of_bins = 5000
min_value = 0
max_value = 3.5

[spacer]
height = 0.15

[E12_Limbs_CTCF_orientation]
file = ${gitHubDirectory}/CTCF/E12_Limbs_Wt_CTCF_colored.bed
title = MACS2_E12_Limbs_Wt_CTCF
height = 0.5
labels = false
line_width = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb

[spacer]
height = 1

[annotations_detailed]
file = ${pathForScripts}/annotations/annotations_detailed.bed
title = annotations detailed
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed

[spacer]
height = 0.5

[Deletion_backgrounds]
file = ${pathForScripts}/annotations/Deletion_backgrounds.bed
title = Deletion_backgrounds
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed
file_type = bed
" > ${ini_file}

pgt --tracks ${ini_file} --region chr2:75122684-75160161 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 36

ini_file="fig2Ctop.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 10000
scalebar_end_position = 97059000
height = 0.3
" > ${ini_file}
for sampleTLA in $samplesTLA; do
  if [[ $sampleTLA = *"CS38" ]]; then
    maxvalue=8500
  else
    maxvalue=5000
  fi
  echo "[$sampleTLA q30 cov]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_q30_cov.bw
title = $sampleTLA mapq30
height = 3
color = black
nans_to_zeros = true
number_of_bins = 25000
min_value = 0
max_value = ${maxvalue}

[spacer]
height = 0.5
" >> ${ini_file}
done
echo "[zoom region]
file = fig2Cbottom.bed
type = vlines
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr10:96980000-97060000 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 20

ini_file="fig2Cbottom.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 100
scalebar_end_position = 97020000
height = 0.3
" > ${ini_file}
for sampleTLA in $samplesTLA; do
  if [[ $sampleTLA = *"CS38" ]]; then
    maxvalue=8500
  else
    maxvalue=5000
  fi
  echo "[$sampleTLA q30 cov]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_q30_cov.bw
title = $sampleTLA mapq30
height = 3
color = black
min_value = 0
max_value = ${maxvalue}

[spacer]
height = 0.5

[$sampleTLA local cov without CATG]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_unaligned_localMapped_onmm10_withoutCATG_cov_UCSC.bdg
title = $sampleTLA local mapped exclude reads containing CATG
height = 3
color = black
min_value = 0
max_value = 600

[spacer]
height = 0.5
" >> ${ini_file}
done
echo "[re]
file = re_fig2Cbottom.bed
display = collapsed
height = 0.5
labels = true
color = black

[zoom region]
file = fig2Cbottom.bed
type = vlines
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr10:97019018-97020046 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 20 \
  --trackLabelFraction 0.3

ini_file="fig3A.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 10000
scalebar_end_position = 75170000
height = 0.3
" > ${ini_file}

maxvalues=('3.5' '2' '1' '2')
i=0
for prot in CTCF RAD21; do
  for genotype in Wt TgN3840; do
    sample=E12_Limbs_${genotype}_${prot}
    if [ $genotype = "TgN3840" ]; then
      file=${pathForChIPGEO}/${sample}_neq2_on_mm10.bw
      macs=${pathForChIPGEO}/${sample}_rep1or2_on_mm10.bed
    else
      file=${pathForChIPGEO}/${sample}_on_mm10.bw
      macs=${pathForChIPGEO}/${sample}.bed
    fi
    echo "[$sample]
file = ${file}
title = ${sample}
height = 3
color = black
nans_to_zeros = true
number_of_bins = 5000
min_value = 0
max_value = ${maxvalues[$i]}

[spacer]
height = 0.15

[MACS2 ${sample}]
file = ${macs}
title = MACS2 ${sample}
color = black
height = 0.2
labels = false
line_width = 1
display = collapsed
file_type = bed

[spacer]
height = 0.15
" >> ${ini_file}
    i=$(($i + 1))
  done
done
echo "[spacer]
height = 0.25

[E12_Limbs_CTCF_orientation]
file = ${gitHubDirectory}/CTCF/E12_Limbs_Wt_CTCF_colored.bed
title = MACS2_E12_Limbs_Wt_CTCF
height = 0.5
labels = false
line_width = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb

[spacer]
height = 1

[annotations_detailed]
file = ${pathForScripts}/annotations/annotations_detailed.bed
title = annotations detailed
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed

[spacer]
height = 0.5


[TgN3840_fosmid]
file = ${pathForScripts}/annotations/TgN3840_fosmid.bed
title = TgN3840_fosmid
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed
file_type = bed

[spacer]
height = 0.5

[Deletion_backgrounds]
file = ${pathForScripts}/annotations/Deletion_backgrounds.bed
title = Deletion_backgrounds
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed
file_type = bed

" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:75105976-75176976 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 36

ini_file="fig3B-E.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 200000
scalebar_end_position = 97800000
height = 0.3
" > ${ini_file}
sample=E12_Limbs_Wt_map_TgN3840
echo "[Hi-C ${sample}]
file = ${pathForHiC}/${sample}/${sample}.40kb.cool
title = $sample 40kb
colormap = ['white', (1, 0.88, 0.66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
depth = 2000000
transform = no
show_masked_bins = true
rasterize = false
min_value = 0
max_value = 0.03

[spacer]
height = 0.25

[E12_Limbs_CTCF_orientation shifted]
file = E12_Limbs_Wt_CTCF_mm10_TgN3840_colored.bed
title = MACS2_E12_Limbs_Wt_CTCF
height = 0.5
labels = false
line_width = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb
" >> ${ini_file}
  for size in 240 480; do
    if [ $size = "240" ]; then
      color=grey
    else
      color=black
    fi
    echo "[spacer]
height = 0.25

[${sample}.${size}kb_domains]
file = ${pathForHiC}/TADs_diff/${sample}.40kb.${size}kb_domains.bed
title = ${sample} 40kb ${size}kb_domains
color = ${color}
border_color = none
height = 0.5
labels = false
line_width = 1
display = interleaved
" >> ${ini_file}
  done

  echo "[spacer]
height = 1

[E12_DFL_Wt_H3K27ac_chr10_shifted_Tg3840]
file = ${pathForChIPGEO}/E12_DFL_Wt_H3K27ac_chr10_shifted_Tg3840.bw
title = E12_DFL_Wt_H3K27ac_chr10_scaled_shifted_Tg3840
height = 3
color = green
nans_to_zeros = true
number_of_bins = 11600
min_value = 0
max_value = 3

[spacer]
height = 1

" >> ${ini_file}
sample=E12_Limbs_TgN3840
for viewpoint in CS38 CS40; do
  echo "[${sample} ${viewpoint}]
file = ${pathFor4CGEO}/${sample}_${viewpoint}.bedGraph.gz
title = ${sample} ${viewpoint}
height = 3
color = red
alpha = 0.5
min_value = 0
max_value = 30
type = line:2
use_middle = true
file_type = bedgraph

[spacer]
height = 0.5
" >> ${ini_file}
done
for viewpoint in CTCF-left CTCF-right; do
  for genotype in Wt TgN3840; do
    sample=E12_Limbs_${genotype}
    echo "[${sample} ${viewpoint}]
file = ${pathFor4CGEO}/${sample}_${viewpoint}.bedGraph.gz
title = ${sample} ${viewpoint}" >> ${ini_file}
    if [ $genotype = Wt ]; then
      echo "height = 3
alpha = 0.8
color = blue" >> ${ini_file}
    else
      echo "overlay_previous = share-y
show_data_range = false
alpha = 0.5
color = red" >> ${ini_file}
    fi
    echo "min_value = 0
max_value = 30
type = line:2
use_middle = true
file_type = bedgraph
" >> ${ini_file}
  done
  echo "[spacer]
height = 0.5
" >> ${ini_file}
done
echo "[4C_viewpoints]
file=${pathFor4CGEO}/../template_4cfile_TgN3840.fa_vp_shifted.bed
title = 4C viewpoints
color = black
height = 0.3
labels = false
line_width = 1
display=collapsed

[spacer]
height = 0.5

[regionsToQuantify4C]
file=${pathForScripts}/annotations/regionsToQuantify4C.bed
title = segments_to_quantify_4C
color = black
height = 0.3
labels = false
line_width = 1
display=collapsed

[spacer]
height = 1

[E12_Limbs_TgN3840_CTCF_neq2_on_TgN3840]
file = ${pathForChIPGEO}/E12_Limbs_TgN3840_CTCF_neq2_on_TgN3840.bw
title = E12_Limbs_TgN3840_CTCF_neq2
height = 3
color = black
nans_to_zeros = true
number_of_bins = 11600
min_value = 0
max_value = 1

[spacer]
height = 0.15

[CTCF orientation]
file = ${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840_colored.bed
color = bed_rgb
border_color = bed_rgb
height = 0.5
display = collapsed
labels = false

[spacer]
height = 1

[E12_Limbs_TgN3840_RAD21_neq2_on_TgN3840]
file = ${pathForChIPGEO}/E12_Limbs_TgN3840_RAD21_neq2_on_TgN3840.bw
title = E12_Limbs_TgN3840_RAD21_neq2
height = 3
color = black
nans_to_zeros = true
number_of_bins = 11600
min_value = 0
max_value = 1

[spacer]
height = 0.15

[E12_Limbs_TgN3840_RAD21_rep1_on_TgN3840]
file = ${pathForChIPGEO}/E12_Limbs_TgN3840_RAD21_rep1_on_TgN3840.bed
title = E12_Limbs_TgN3840_RAD21_rep1_on_TgN3840
color = black
height = 0.2
labels = false
line_width = 1
display = collapsed
file_type = bed

[E12_Limbs_TgN3840_RAD21_rep2_on_TgN3840]
file = ${pathForChIPGEO}/E12_Limbs_TgN3840_RAD21_rep2_on_TgN3840.bed
title = E12_Limbs_TgN3840_RAD21_rep2_on_TgN3840
color = black
height = 0.2
labels = false
line_width = 1
display = collapsed
overlay_previous = yes
file_type = bed

[spacer]
height = 0.5

[Genes]
# file = custom_mm10_TgN3840_protein_coding.gtf
file = gencode.vM23.protcod.mm10_TgN3840.gtf
title = genes on mm10_TgN3840_protein_coding
color = black
height = 0.5
labels = false
line_width = 0.5
display = collapsed
merge_transcripts = true

[integration]
file = integration_borders_TgN3840.bed
type = vlines
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr10:95480001-97880000 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 34 \
  --trackLabelFraction 0.3


ini_file="fig4A-C.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 200000
scalebar_end_position = 97800000
height = 0.3
" > ${ini_file}

for genotype in Wt TgN3840; do
  sample=E12_Limbs_${genotype}_map_TgN3840
  echo "[Hi-C ${sample}]
file = ${pathForHiC}/${sample}/${sample}.40kb.cool
title = $sample 40kb
colormap = ['white', (1, 0.88, 0.66), (1, 0.74, 0.25), (1, 0.5, 0), (1, 0.19, 0), (0.74, 0, 0), (0.35, 0, 0)]
depth = 2000000
transform = no
show_masked_bins = true
rasterize = false
min_value = 0
max_value = 0.03

[spacer]
height = 0.25
" >> ${ini_file}
  if [ $genotype = Wt ]; then
    file=E12_Limbs_Wt_CTCF_mm10_TgN3840_colored.bed
    title=E12_Limbs_Wt
  else
    file=${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840_colored.bed
    title=E12_Limbs_TgN3840_CTCF_rep1or2
  fi
    echo "[${title}]
file = ${file}
title = ${title}
height = 0.5
labels = false
line_width = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb

[spacer]
height = 0.25" >> ${ini_file}
  for size in 240 480; do
    if [ $size = "240" ]; then
      color=grey
    else
      color=black
    fi
    echo "[spacer]
height = 0.25

[${sample}.${size}kb_domains]
file = ${pathForHiC}/TADs_diff/${sample}.40kb.${size}kb_domains.bed
title = ${sample} 40kb ${size}kb_domains
color = ${color}
border_color = none
height = 0.5
labels = false
line_width = 1
display = interleaved
" >> ${ini_file}
  done

  echo "[spacer]
height = 0.25
" >> ${ini_file}
  size=240
    echo "
[${sample}.${size}kb_tad_score]
file = ${pathForHiC}/TADs_diff/${sample}.40kb.${size}kb_tad_score.bm
title = ${sample} 40kb all tad_score
color = black
use_middle = true
min_value = -2
max_value = 2
type = line:1
height = 3
file_type = bedgraph
" >> ${ini_file}
  for size in 480; do
    echo "[${sample}.${size}kb_tad_score]
file = ${pathForHiC}/TADs_diff/${sample}.40kb.${size}kb_tad_score.bm
color = grey
use_middle = true
type = line
overlay_previous = share-y
file_type = bedgraph" >> ${ini_file}
  done

  echo "[hline]
file_type = hlines
y_values = 0
overlay_previous = share-y
line_style = dashed

[spacer]
height = 0.25
" >> ${ini_file}
done
echo "[Difference 40kb]
file = ${pathForHiC}/TADs_diff/E12_Limbs_TgN3840MinusWt_map_TgN3840.40kb.cool
title = E12_Limbs_TgN3840MinusWt_map_TgN3840 40kb
colormap = bwr
depth = 2000000
transform = no
show_masked_bins = true
rasterize = false
min_value = -1.5e-7
max_value = 1.5e-7

[Quantification_region]
file = interSubTAD.bedpe
file_type = links
links_type = loops
line_width = 2
color = black
line_style = dashed
overlay_previous = share-y

[spacer]
height = 1

[Genes]
# file = custom_mm10_TgN3840_protein_coding.gtf
file = gencode.vM23.protcod.mm10_TgN3840.gtf
title = genes on mm10_TgN3840_protein_coding
color = black
height = 0.5
labels = false
line_width = 0.5
display = collapsed
merge_transcripts = true

[integration]
file = integration_borders_TgN3840.bed
type = vlines
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr10:95480001-97880000 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 34 \
  --trackLabelFraction 0.3

ini_file="figS1A.ini"
sampleTLA=TgN3840_vector
echo "[$sampleTLA q30 cov]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_q30_cov.bw
title = $sampleTLA mapq30
height = 3
color = black
nans_to_zeros = true
number_of_bins = 1955
summary_method = max
min_value = 0
max_value = 5000

[spacer]
height = 0.5

[x-axis]
" > ${ini_file}

for chr in {1..19} X Y; do
  size=$(cat $pathForChromSizes | awk -v chr=$chr '$1=="chr"chr{print $2}')
  # We want that 200Mb = 89cm
  pdf_size=$(echo $size | awk '{print $1 * 89 / 200000000}')
  echo "chr${chr}:0-${size}"
  pgt --tracks ${ini_file} --region chr${chr}:0-${size} \
    -o ${ini_file/.ini/_chr${chr}_from0.pdf} --dpi 500 --plotWidth ${pdf_size}
done

ini_file="figS1B.ini"
sampleTLA=TgN3840_vector
echo "[$sampleTLA onTgN3840_fosmid_TH_TT_HH]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_onTgN3840_fosmid_TH_TT_HH.bw
title = $sampleTLA onTgN3840_fosmid_TH_TT_HH
height = 3
color = black
nans_to_zeros = true
number_of_bins = 1955
summary_method = max
min_value = 0
max_value = 1000

[spacer]
height = 0.5

[x-axis]
" > ${ini_file}
pgt --tracks ${ini_file} --region TgN3840_fosmid_TH_TT_HH:0-44996 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 20 \
  --trackLabelFraction 0.3

ini_file="figS1C.ini"
sampleTLA=TgN3840_vector
echo "[$sampleTLA onTgN3840_fosmid_TH_TT_HH]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_onTgN3840_fosmid_TH_TT_HH.bw
title = $sampleTLA onTgN3840_fosmid_TH_TT_HH
height = 3
color = black
nans_to_zeros = true
number_of_bins = 1955
summary_method = max
min_value = 0
max_value = 1000

[spacer]
height = 0.5

[$sampleTLA onTgN3840_fosmid_TH_TT_HH]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_onTgN3840_fosmid_TH_TT_HH_summary.bed
title = $sampleTLA onTgN3840_fosmid_TH_TT_HH as bed
height = 10
color = black
labels = true

[spacer]
height = 0.5

[vlines]
file = fosmid_boundaries_on_TH_TT_HH.bed
type = vlines

[x-axis]
" > ${ini_file}
pgt --tracks ${ini_file} --BED fosmid_boundaries_on_TH_TT_HH_1kb.bed \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 20 \
  --trackLabelFraction 0.3

ini_file="figS2B.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 5000
center = 75133816
height = 0.3

[spacer]
height = 0.25
" > ${ini_file}
sample=E12_Limbs_TgN3840_input
for ws in 1 2; do
  for prof in ratio meanRatio; do
    echo "[$ws $prof]
file = ${pathForFREECResults}/${sample}_${ws}kb_${prof}.bedgraph
title = $prof $ws kb
height = 3
color = black
nans_to_zeros = true
min_value = 0
max_value = 6
" >> ${ini_file}
    if [ $prof = "meanRatio" ]; then
      echo "[CN_expected_TLA]
file = ${pathForScripts}/annotations/expected_CN_TLA_ignore27bpCS39.bedgraph
height = 3
type = line
color = #be1e2d
overlay_previous = share-y
show_data_range = false

[CN_expected_TLA+1copy]
file = ${pathForScripts}/annotations/expected_CN_TLA+1copy_ignore27bpCS39.bedgraph
height = 3
type = line
color = #fdba13
overlay_previous = share-y
show_data_range = false
" >> ${ini_file}
    fi
    echo "[spacer]
height = 0.15
" >> ${ini_file}
  done
done

echo "[CTCF orientation]
file = ${gitHubDirectory}/CTCF/E12_Limbs_TgN3840_CTCF_rep1or2_on_TgN3840_colored.bed
color = bed_rgb
border_color = bed_rgb
height = 0.5
display = collapsed
labels = false
[spacer]
height = 1

[annotations_detailed]
file = ${pathForScripts}/annotations/annotations_detailed.bed
title = annotations detailed
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed

[spacer]
height = 0.5


[TgN3840_fosmid]
file = ${pathForScripts}/annotations/TgN3840_fosmid.bed
title = TgN3840_fosmid
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed
file_type = bed

[spacer]
height = 0.5

[Deletion_backgrounds]
file = ${pathForScripts}/annotations/Deletion_backgrounds.bed
title = Deletion_backgrounds
color = black
height = 0.5
labels = true
line_width = 1
display = collapsed
file_type = bed

" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:75105976-75176976 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 21.852


ini_file="figS3A.ini"
echo "[x-axis]

[scalebar]
file_type = scalebar
size = 100
scalebar_end_position = 97020000
height = 0.3
" > ${ini_file}
for sampleTLA in $samplesTLA; do
  if [[ $sampleTLA = *"CS38" ]]; then
    maxvalue=8500
    maxvalueLocal=600
  else
    maxvalue=5000
    maxvalueLocal=5000
  fi
  echo "[$sampleTLA q30 cov]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_q30_cov.bw
title = $sampleTLA mapq30
height = 3
color = black
min_value = 0
max_value = ${maxvalue}

[spacer]
height = 0.5

[$sampleTLA local cov]
file = ${pathForTLA}/${sampleTLA}/${sampleTLA}_unaligned_CATG_split_nomm10noTgN3840_fosmid_localMapped_onmm10.bw
title = $sampleTLA CATG split no mm10 no fosmid local mapped
height = 3
color = black
min_value = 0
max_value = ${maxvalueLocal}

[spacer]
height = 0.5
" >> ${ini_file}
done
echo "[re]
file = re_fig2Cbottom.bed
display = collapsed
height = 0.5
labels = true
color = black

[integration]
file = integration_borders.bed
type = vlines
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr10:97019018-97020046 \
  -o ${ini_file/ini/pdf} --dpi 500 --plotWidth 20 \
  --trackLabelFraction 0.3

### STATISTICS
# On 4C:
Rscript ${pathForScripts}/4C_quantification.R \
  ${pathFor4CGEO}/segToFrag_E12_Limbs_TgN3840_CTCF-left.bw \
  ${pathFor4CGEO}/segToFrag_E12_Limbs_Wt_CTCF-left.bw \
  ${pathForScripts}/annotations/regionsToQuantify4C.bed \
  welc_TAD_left_of_int_mm10_613K_coordinates \
  welc_TAD_right_of_int_mm10_613K_coordinates 2> quantif.log
# left of integration: +11%
# right of integration: -47%
Rscript ${pathForScripts}/4C_quantification.R \
  ${pathFor4CGEO}/segToFrag_E12_Limbs_TgN3840_CTCF-right.bw \
  ${pathFor4CGEO}/segToFrag_E12_Limbs_Wt_CTCF-right.bw \
  ${pathForScripts}/annotations/regionsToQuantify4C.bed \
  welc_TAD_left_of_int_mm10_613K_coordinates \
  welc_TAD_right_of_int_mm10_613K_coordinates 2>> quantif.log
# left of integration: -32%%
# right of integration: +64%%

# On Hi-C
sample1=E12_Limbs_TgN3840_map_TgN3840
cool1=${pathForHiC}/${sample1}/${sample1}.40kb.cool
sample2=E12_Limbs_Wt_map_TgN3840
cool2=${pathForHiC}/${sample2}/${sample2}.40kb.cool
python ${pathForScripts}/quantifyCooler.py \
  --bedpe interSubTAD.bedpe --cool1 ${cool1}\
  --cool2 ${cool2}

# Region  Mean_first_file Mean_second_file        ratio   pval_Mann_Whitney_U_test        pval_Wilcoxon_signed_rank_test
# interSubTAD     0.0075496819956354265   0.012972878417607676    0.5819588955207106      1.8909929157763953e-32  1.5259922107013315e-25

# Boundaries
Rscript ${pathForScripts}/getBoundariesTable.R ${pathForHiC} ${pathForScripts}/annotations/annotated_boundaries.bed boundaries_table.txt

# Figure 2e:
cd ${pathMinION}

bash ${gitHubDirectory}/scripts/minION/splitBigFasta.sh targetted_reads.fa

ln -s ${gitHubDirectory}/minION_nCATS/TgN3840_TLAderived.fa .

wget "jura.wi.mit.edu/page/papers/Hughes_et_al_2005/tables/dot_plot.pl" -O ${gitHubDirectory}/scripts/minION/dot_plot.pl

# Generate the dotplot, this is slow:
Rscript ${gitHubDirectory}/scripts/minION/dotplot_selected_reads.R

mv fig2E.pdf ${path}
