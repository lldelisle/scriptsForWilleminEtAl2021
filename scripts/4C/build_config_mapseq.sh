echo "[Job]"
newDesc=`echo $1 | awk -F '/' '{print "mappingAfter"$(NF-2)}'`
echo "description='$newDesc'"
echo "assembly_id='mm10'"
echo ""
echo "[Options]"
echo "input_type=0 #align on the genome. Set to 1 for exonome and 2 to transcriptome"
echo "compute_densities=True"
echo "merge_strands=0 "
echo "discard_pcr_duplicates=False"
echo "create_gdv_project=False #always leave it to False"
echo ""
echo "[Groups]"
i=1
for f in "$@"; do
  echo "[[${i}]]"
  groupName=`basename $f | awk '{gsub(".fastq$","",$1);gsub("_filtered$","",$1);print $1}'`
  echo "name='${groupName}'"
  echo "control=False"
  let i++
done
echo ""
echo "[Runs]"
i=1
for f in "$@"; do
  echo "[[${i}]]"
  echo "url='$f'"
  echo "group_id=$i"
  let i++
done
