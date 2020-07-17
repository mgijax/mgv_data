#
# GffSorter.py
#
# Sorts the contents of a GFF3 file so that:
#
# 1. All features for a chromosome are grouped (no intervening features from other chromosomes)
# 2. All features for a gene model are grouped (no intervening features from other genes)
# 3. The order of gene models for a chromosome is by start position of the top-level (gene) feature.
# 4. Within a gene model, features are ordered so there are no forward references.
# 5. There is a "###" terminator after each gene group
#

import sys
try:
    from .gff3lite import Gff3Parser, formatLine
except:
    from gff3lite import Gff3Parser, formatLine

def sorted (inGff) :
     return Gff3Parser(inGff, returnHeader = True, returnGroups=True).sortIterate()

def sort (inGff, outGff) :
    sortedStream = sorted(inGff)
    header = next(sortedStream)
    for l in header:
        outGff.write(l)
    for grp in sortedStream:
        for f in grp:
            outGff.write(formatLine(f))
        outGff.write("###\n")

if __name__ == "__main__":
    sort(sys.stdin, sys.stdout)
