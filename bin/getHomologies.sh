#!/bin/bash
#
# Computes homology clusters (connected components) from Alliance data, and assigns 
# cluster ids to gene ids.
# Downloads three set of homology pairs from the Alliance (mouse-human, mouse-rat,
# and human-rat), builds a graph, and enumerates its connected components. 
#

source config.sh
source utils.sh

taxonIds=()
until [ -z "$1" ]  # Until all parameters used up . . .
do
    case "$1" in
    -x)
        # set the organism path name, eg mus_musculus_dba2j
        shift
        taxonIds+=("$1")
        ;;
    *)
        # anything else is an error
        echo "Unrecognized option:" $1
        usage
    esac
    shift
done

nt=${#taxonIds[@]}

function getAlliancePairs {
  ga="$1"
  gb="$2"
  fname="$3"
  logit "Get alliance pair: $1 $2"
  if [ ! -f "$fname" ]; then
      curl -X GET "https://www.alliancegenome.org/api/homologs/${ga}/${gb}?stringencyFilter=stringent&rows=200000&start=1" -H "accept: application/json" > "${fname}"
  else
      logit "Warning: Skipping download. File exists: ${fname}"
  fi
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
mkdir -p ${DDIR}
logit "${PYTHON} getHomologies.py ${files[@]}"
${PYTHON} getHomologies.py ${files[@]}
#
#rm ${files[@]}
