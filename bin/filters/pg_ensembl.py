#
# Munges a feature from an Ensembl gff3 for import to MGV.
# 1. for genes, sets col 3 type to "protein_coding_gene" as appropriate.
#
def feature (f) :
  attrs = f[8]
  if f[2] == "gene" and attrs.get('biotype', None) == 'protein_coding':
    f[2] = 'protein_coding_gene'
  return f
