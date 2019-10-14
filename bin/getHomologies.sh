#!/bin/bash
#
# Computes mouse-human-rat homology clusters from Alliance data, and assigns 
# cluster ids to gene ids.
# Downloads three set of homology pairs from the Alliance (mouse-human, mouse-rat,
# and human-rat), builds a graph, and enumerates its connected components. 
#

source config.sh
source utils.sh

taxonIds=($*)
nt=${#taxonIds[@]}

function getAlliancePairs {
  ga="$1"
  gb="$2"
  fname="$3"
  curl -X GET "https://www.alliancegenome.org/api/homologs/${ga}/${gb}?stringencyFilter=stringent&rows=200000&start=1" -H "accept: application/json" > "${fname}"
}

for (( c=0; c<${nt}; c++ ))
do  
  for (( d=${c}; d<${nt}; d++ ))
  do  
     if [ ${c} != ${d} ] ; then
       t1=${taxonIds[${c}]}
       t2=${taxonIds[${d}]}
       file="${DDIR}/${t1}-${t2}.json "
       getAlliancePairs $t1 $t2 $file
       files+="${file}"
     fi
  done
done
#
python getHomologies.py ${files[@]} > ${DDIR}/homologies.txt
#
rm ${files[@]}
