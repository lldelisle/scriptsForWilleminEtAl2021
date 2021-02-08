#!/bin/bash

inputFasta=$1


if [[ $inputFasta =~ \.gz$ ]]; then
  cmd="zcat"
else
  cmd="cat"
fi

eval "$cmd $inputFasta | csplit -z -s -b %09d -f "${inputFasta}_" - '/>/' '{*}' &"
wait

mkdir -p ${inputFasta}CONTIGS
for f in ${inputFasta}_*;
  do 
  n=`head -n 1 $f | awk '{print substr($1,2)}'`
  mv $f ${inputFasta}CONTIGS/${n}.fa
done
