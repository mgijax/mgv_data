#
# pg_MGI.py
#
# Munges a record from the MGI gff3 file for MGV.
# For genes and pseudogenes (top level features):
# 1. Moves the so_term_name attribute into column 3
# 2. Moves the curie attribute to the cID attribute.
# 3. Moves transcript and CDS ids into ID attrs
# All other features are unaffected.
#

import sys
def header (header):
    attrs = {
    'genome-build' : 1,
    'genome-version' : 1,
    'genome-date' : 1,
    'genome-build-accession' : 1,
    'genebuild-last-updated' : 1,
    }
    h2 = []
    for l in header:
      if l.startswith('#!'):
        n = l.split()[0][2:]
        v = attrs.pop(n, 0)
        if v:
          h2.append(l)
    return h2

i2i = {}
def feature(f):
    attrs = f[8]
    if f[2] in ["gene","pseudogene"]:
        #
        f[2] = attrs['so_term_name']
        del attrs['so_term_name']
        #
        attrs['cID'] = attrs['curie']
        del attrs['curie']
    elif f[2] == "CDS":
        pid = attrs.get('protein_id', attrs['ID'])
        attrs['ID'] = pid
    elif f[2] != "exon":
        tid = attrs.get('transcript_id', attrs['ID'])
        i2i[attrs['ID']] = tid
        attrs['ID'] = tid
    if 'Parent' in attrs:
        attrs['Parent'] = map(lambda p: i2i.get(p, p), attrs['Parent'])
    return f
