#
# pg_MGI.py
#
# Munges a record from the MGI gff3 file for MGV.
# For genes and pseudogenes (top level features):
# 1. Moves the so_term_name attribute into column 3
# 2. Moves the curie attribute to the cID attribute.
# All other features are unaffected.
#

def feature(f):
    if f[2] in ["gene","pseudogene"]:
        attrs = f[8]
        #
        f[2] = attrs['so_term_name']
        del attrs['so_term_name']
        #
        attrs['cID'] = attrs['curie']
        del attrs['curie']
    return f
