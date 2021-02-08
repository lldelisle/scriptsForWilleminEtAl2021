wd=$1
genome=$2
assembly=$3
pathWithMutantGenome=$4
templatePrimerFile=$5
resultDir=${wd}/res_files_mapping_${genome}
echo "[Job]" > config_4cseq_${genome}.txt
newDesc="4Cafter$(basename $wd)_on${assembly}"
echo "description='$newDesc'" >> config_4cseq_${genome}.txt
echo "assembly_id='$assembly'" >> config_4cseq_${genome}.txt
echo "" >> config_4cseq_${genome}.txt
echo "[Options]" >> config_4cseq_${genome}.txt
echo "#merge_strands=1 !! only define if no url_wig provided in [Groups]. Can be set to 0 otherwise" >> config_4cseq_${genome}.txt
echo "primersfile='$wd/4c_primerFile_${genome}.fa'" >> config_4cseq_${genome}.txt
echo "norm_reg=5000000" >> config_4cseq_${genome}.txt
echo "" >> config_4cseq_${genome}.txt
echo "[Groups]" > config_4cseq_${genome}_groups.txt
echo "[Runs]" > config_4cseq_${genome}_runs.txt
if [ -e $wd/4c_primerFile_${genome}.fa ]; then
  rm $wd/4c_primerFile_${genome}.fa
fi
j=1
lastj=0
lastgroupname=""
chunk=1
maxChunk=7
groupInChunk=0
for f in ${resultDir}*/*.bam ; do
  sample=`basename $f | awk '{gsub(".bam$","",$1);gsub("_filtered$","",$1);print $1}'`
  # If there are replicates they will be merged into one group:
  groupName=`echo $sample | awk '{split($1, r, "_r"); print r[1]}'`
  # If it is a second/third... replicate we do not put it into the group list
  if [ "$groupName" != "$lastgroupname" ]; then
    if [ $groupInChunk -eq $maxChunk ]; then
      # We need to change chunk:
      # First we write down the current config_4cseq
      cat config_4cseq_${genome}.txt config_4cseq_${genome}_groups.txt config_4cseq_${genome}_runs.txt > config_4cseq_${genome}_split${chunk}.txt
      # Update the description
      sed -i "s/${newDesc}/${newDesc}_${chunk}/" config_4cseq_${genome}_split${chunk}.txt
      # We prepare the new files
      echo "[Groups]" > config_4cseq_${genome}_groups.txt
      echo "[Runs]" > config_4cseq_${genome}_runs.txt
      # We update the variables
      let chunk++
      groupInChunk=0
    fi
    lastgroupname=$groupName
    lastj=$j
    let groupInChunk++
    n=`ls ${resultDir}*/${groupName}*.bam | wc -l`
    if [ $n -gt 1 ]; then
      groupName="${groupName}_neq${n}"
    fi
    echo "[[${j}]]" >> config_4cseq_${genome}_groups.txt
    echo "name='${groupName}'" >> config_4cseq_${genome}_groups.txt
    echo "library_file_type_id=''" >> config_4cseq_${genome}_groups.txt
    echo "library_id=''" >> config_4cseq_${genome}_groups.txt
    echo "library_file_url='${pathWithMutantGenome}/${genome}/library_${genome}_NlaDpn_30bps_segmentInfos.bed' #set the path the the lib" >> config_4cseq_${genome}_groups.txt
    echo "library_param_file=''" >> config_4cseq_${genome}_groups.txt
    echo "window_size=11" >> config_4cseq_${genome}_groups.txt
    echo "run_domainogram=False" >> config_4cseq_${genome}_groups.txt
    echo "before_profile_correction=False" >> config_4cseq_${genome}_groups.txt
    echo "" >> config_4cseq_${genome}_groups.txt
    # We add the group to the primer file.
    # We assume that the name of the viewpoint is part of the groupName separated by _ and that the viewpoint is part of the template.
    cat ${templatePrimerFile} | awk -F "|" -v gn=${groupName} -v OFS="|" '
      BEGIN{
        split(gn, gns, "_")
        for (i in gns) keywords[gns[i]] = ""
        written = 0
      }
      written == 0 && (substr($1, 2) in keywords){
        $1 = ">"gn
        print $0
        getline
        print $0
        written = 1
      }
      END{
        if(written == 0){
          print ">"gn"|notFound"
          print "ATCG"
        }
      }' >> $wd/4c_primerFile_${genome}.fa
      if [ `grep notFound $wd/4c_primerFile_${genome}.fa | wc -l` -ne 0 ]; then
        echo "There is a error in a viewpoint"
        grep notFound $wd/4c_primerFile_${genome}.fa
        exit 1
      fi
  fi
  # We also write the Runs part:
  echo "[[${j}]]" >> config_4cseq_${genome}_runs.txt
  echo "url='${f}'" >> config_4cseq_${genome}_runs.txt
  wig=`ls ${resultDir}*/${sample}_merged.sql`
  echo "url_wig='${wig}'" >> config_4cseq_${genome}_runs.txt
  echo "group_id=${lastj}" >> config_4cseq_${genome}_runs.txt
  let j++
done
# We now write the last config file:
cat config_4cseq_${genome}.txt config_4cseq_${genome}_groups.txt config_4cseq_${genome}_runs.txt > config_4cseq_${genome}_split${chunk}.txt
# Update the description
sed -i "s/${newDesc}/${newDesc}_${chunk}/" config_4cseq_${genome}_split${chunk}.txt
# And delete the temporary files
rm config_4cseq_${genome}.txt config_4cseq_${genome}_groups.txt config_4cseq_${genome}_runs.txt

