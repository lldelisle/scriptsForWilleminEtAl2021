brFiles=$1
primerFile=$2
pathWithBRFiles=$3
pathForScripts=$4

# Extract from the 4cfile the position of the viewpoint
# and put it into a bed
cat $primerFile | awk -F "|" -v OFS="\t" 'NR%2==1{split($3,a,":|-");print a[1],a[2],a[3],$1,1,"+"}' > ${primerFile}_vp.bed

# We shift the bed 
if [[ "$brFiles" = *".R" ]]; then
  Rscript ${pathWithBRFiles}${brFiles} ${primerFile}_vp.bed 1 2 3 6 ${primerFile}_vp_shifted.bed
else
  brFilesSpace=$(echo $brFiles | tr "," " ")
  for f in $brFilesSpace; do
    ln -s ${pathWithBRFiles}$f .
  done
  bash $pathForScripts/shiftBedWithMultipleBR.sh $brFiles ${primerFile}_vp.bed ${primerFile}_vp_shifted.bed 1
fi

mv $primerFile ${primerFile}_ori
# The primer file is formed of groups of 2 lines.
# The first line is description
# The second line is sequence
# We only modify the first lines
h=1
while read line; do
  if [ $h = 1 ]; then
    # We look into the shifted bed which line corresponds to the viewpoint
    # Then we correct the line with the new viewpoint position
    # We also correct the exclude part using the same distance on each side
    cat ${primerFile}_vp_shifted.bed | awk -v l=$line '
      BEGIN{
        split(l,ls,"|")
        found=0
      }
      $4==ls[1]&&found==0{
        oldVP=ls[3]
        current_chr=$1
        current_start=$2
        current_end=$3
        found=1
      }
      $4==ls[1]&&found==1{
        if ($1 == current_chr){
          if ($2 < current_start){
            current_start = $2
          }
          if ($3 > current_end){
            current_end = $3
          }
        } else {
          current_chr="INCONSISTENT"
          current_start = ""
          current_end = ""
        }
      }
      END{
        ls[3]=current_chr":"current_start"-"current_end
        for(i=5;i<=length(ls);i++){
          if(gsub("Exclude=","",ls[i])>0){
            split(ls[i],lse,":|-")
            split(oldVP,lsvp,":|-")
            if($6=="+"){
              ls[i]="Exclude="current_chr":"current_start-(lsvp[2]-lse[2])"-"current_end+(lse[3]-lsvp[3])
            }else{
              ls[i]="Exclude="current_chr":"current_start-(lse[3]-lsvp[3])"-"current_end+(lsvp[2]-lse[2])
            }
          }
        }
        printf("%s",ls[1])
        for(i=2;i<=length(ls);i++){
          printf("|%s",ls[i])
        }
        print ""
      }'  >> $primerFile
    h=0
  else
    echo $line >> $primerFile
    h=1
  fi
done < ${primerFile}_ori
