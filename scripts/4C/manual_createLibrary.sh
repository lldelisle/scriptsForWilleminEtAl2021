#!/bin/bash

# create library:
# ./createLibrary -i mm9 -m CATG -s GATC -l 30 -n NlaDpnLibrary_30bps
echo "author: Marion Leleu"
echo `date`

scriptsDirectory="./"
bedToolsDirectory=""
echo "will call script ${scriptsDirectory}getRestEnzymeOccAndSeq.pl"

arguments=$*
echo "arguments:"${arguments}

while getopts i:g:m:s:l:n:r:w: name  
  do
    case $name in
        i)
            genomeFile="$OPTARG"
            ;;
    g)
        genomeName_flag=1
        genomeName="$OPTARG"
        ;;
        m)
            mainRestrictionSite="$OPTARG"
            ;;
        s)
            secondaryRestrictionSite="$OPTARG"
            ;;
        l)
            segmentLength="$OPTARG"
            ;;
        n)
            libraryName="$OPTARG"
            ;;
        r)
            repeatsFile_flag=1
            repeatsFile="$OPTARG"
            ;;
        w)
            wd="$OPTARG"
            ;;
        ?)
            printf "Usage: %s: [-i genome file (fasta, incl path)] [-g genome name (required if no genome file)] [-m primary restriction site] [-s secondary restriction site] [-l segment length] [-n library name (default=myLibrary] [-r repeat masker bed file (if not provided then use the genome name to find it in existing ones)] [-w working dir]\n" $0
            exit 2
    esac
done

if [[ -z ${wd} ]]
then
    wd=`pwd`
    wd=${wd}"/"
fi
echo "wd="${wd}

if [[ -z ${genomeFile} && -z ${genomeName_flag} ]]
then
    echo "No genome file or genome name provided. At least one is required"
    echo "Creation of the library has been cancelled"
        exit 2

fi

if [[ -z ${genomeFile} && ! -z ${genomeName_flag} ]]
then
    genomeIndex=`awk -v genome=${genomeName} '{if($1 ~ genome){print $2".fa"}}' ${scriptsDirectory}genomesCorrespondanceTable.txt`
    if [[ ! -e ${genomeIndex} ]]
    then
            echo "No genome index found for "${genomeName}
            echo "Creation of the library has been cancelled"
        exit 2
    else
        genomeFile=${genomeIndex}
    fi
fi


if [[ `echo ${mainRestrictionSite} | awk '{gsub(/[ATGC]/,"",$0);if(length($0)>0){print "not OK"}}'` == "not OK" ]]
then
    echo "Invalid primary restriction site"
    echo "Creation of the library has been cancelled"
    exit;
fi

if [[ ${secondaryRestrictionSite} && `echo ${secondaryRestrictionSite} | awk '{gsub(/[ATGC]/,"",$0);if(length($0)>0){print "not OK"}}'` == "not OK" ]]
then
        echo "Invalid secondary restriction site"
        echo "Creation of the library has been cancelled"
        exit;
fi

if [[ -z ${segmentLength} ]]
then
    segmentLength=30
fi

if [[ -z ${libraryName} ]]
then
    libraryName="myLibrary"
fi

segmentFile="${wd}${libraryName}.txt"
fragmentFile="${wd}${libraryName}.frag"

cmdPerl="${scriptsDirectory}getRestEnzymeOccAndSeq.pl -i ${genomeFile} -m ${mainRestrictionSite} -s ${secondaryRestrictionSite} -l ${segmentLength} -o ${segmentFile} -f ${fragmentFile}"
echo "cmdPerl="${cmdPerl}
perl ${cmdPerl}

segmentBedFile=${wd}${libraryName}".bed"
fragmentBedFile=${wd}${libraryName}"_frag.bed"
segInfoBedFile=${wd}${libraryName}"_segmentInfos.bed"
segInfoBedFileWithRepeatsCoverage=${wd}${libraryName}"_segmentInfos_repeatsCoverage.bed"
echo -e  "\nPrepare bed files"
awk '{if(NR>1){chr=$2; gsub("Mt","M",chr);print chr"\t"$3"\t"$4"\tfrag"$1;}}' ${fragmentFile} > ${fragmentBedFile}
awk '{if(NR>1){chr=$2; gsub("Mt","M",chr);if($9<0){$9=0};print chr"\t"$6"\t"$7"\tfrag"$1"_startSeq"; print chr"\t"$9"\t"$10"\tfrag"$1"_endSeq";}}' ${fragmentFile} > ${segmentBedFile}
awk '{if(NR>1){fragmentInfo=$1"|"$2":"($3+1)"-"$4"|indexOfSecondRestSiteOcc="$11"|status="$NF"|length="$4-$3"|0|0|0|0";if($9<0){$9=0};print $2"\t"$6"\t"$7"\ttype=startSegment|"fragmentInfo;print $2"\t"$9"\t"$10"\ttype=endSegment|"fragmentInfo;}}' ${fragmentFile} > ${segInfoBedFile}


if [[ -z ${repeatsFile_flag} ]]
then
    repeatsFile="${repeatsPath}/${genomeName}/${genomeName}_rmsk.bed"
fi

if [[ -f ${repeatsFile} ]]
then
    echo -e "\nCalculate repeats coverage with file ${repeatsFile}"

    ${bedToolsDirectory}coverageBed -a ${segInfoBedFile} -b ${repeatsFile} | awk '{n=split($4,a,"|"); infos=""; for(i=1;i<=n-4;i++){infos=infos""a[i]"|"} print $1"\t"$2"\t"$3"\t"infos""$5"|"$6"|"$7"|"$8;}' | sort -k1,1 -k2,3n > ${segInfoBedFileWithRepeatsCoverage}

    mv ${segInfoBedFileWithRepeatsCoverage} ${segInfoBedFile}
else
    echo -e "\n Did not calculate repeats coverage as not rmsk file found nor provided"
fi

grep "FragIsValid" ${segInfoBedFile} > ${segInfoBedFile}".tmp"
mv ${segInfoBedFile}".tmp" ${segInfoBedFile}


echo "library files created!"
echo "main file: ${segInfoBedFile}"
