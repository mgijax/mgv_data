#!/bin/bash
downloadsdir=${1}
outputdir=${2}
chunksize=4000000
# Part 1. Generate gff3 files from genomes in MouseMine
python getGenomeFromMouseMine.py -o ${downloadsdir}
# Part 2. Import genomes - chunk the files.
for  i in $( ls ${downloadsdir} ); do
  python importGenome.py -d ${outputdir} -k ${chunksize} < ${downloadsdir}/$i
done
