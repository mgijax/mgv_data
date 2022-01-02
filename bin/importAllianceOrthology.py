#
# importAllianceOrthology.py
#
# Reads from stdin. Writes multiple output files to ${ODIR}/homologies
#
# Extracts the essential bits from the Alliance strict orthology set (the IDs).
# Writes one file per taxon id containing the associations for that taxon.
#
# Input is tab delimited with these columns:
# Gene1ID
# Gene1Symbol
# Gene1SpeciesTaxonID
# Gene1SpeciesName
# Gene2ID
# Gene2Symbol
# Gene2SpeciesTaxonID
# Gene2SpeciesName
# Algorithms
# AlgorithmsMatch
# OutOfAlgorithms
# IsBestScore
# IsBestRevScore
#
# Output includes (Gene1ID, Gene1SpeciesTaxonID, Gene2ID, Gene2SpeciesTaxonID)
# ("NCBITaxon:" is stripped off fields 2 and 4)
#
# Example line:
#     HGNC:18487      ST13P4  NCBITaxon:9606  Homo sapiens    FB:FBgn0260484  HIP     NCBITaxon:7227  Drosophila melanogaster Hieranoid|OrthoInspector|OrthoFinder|InParanoid|PANTHER     5       10      Yes     No
# This would add
#     HGNC:18487      9606  FB:FBgn0260484  7227
# to the output file 9606.tsv
#

import sys
import os

odir = os.path.join(os.environ["ODIR"],"homologies")
tx2file = {}
first = True
for line in sys.stdin:
    # skip comment lines
    if line.startswith("#"):
        continue
    # skip header line
    if first:
        first = False
        continue
    #
    fields = line[:-1].split("\t")
    id1, tx1, id2, tx2 = fields[0], fields[2], fields[4], fields[6]
    tx1 = tx1.replace('NCBITaxon:', '')
    tx2 = tx2.replace('NCBITaxon:', '')

    oline = "\t".join([id1, tx1, id2, tx2]) + "\n"
    ofile = os.path.join(odir, tx1) + ".tsv"
    if not tx1 in tx2file:
       tx2file[tx1] = open(ofile, 'w')
    ofd = tx2file[tx1]
    ofd.write(oline)

