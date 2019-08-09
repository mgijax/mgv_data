#!/bin/bash

source utils.sh

while read -u 10 p; do
  IFS='	' read -r -a array <<< "${p}"
  srelease="${array[0]}" 
  if [[ ${srelease} != "#"* ]] ; then
    staxon="${array[1]}" 
    spath="${array[2]}" 
    sname="${array[3]}" 
    ./getGenome.sh -r ${srelease} -g ${spath} -n ${sname} -x ${staxon} --mgi-models $*
    checkExit
  fi
done 10<${GENOMES_FILE}
