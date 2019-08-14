#
# 1. Removes "gene:" "transcript:" and "cds:" prefixes from IDs and Parents.
#
def stripPrefix(eid) :
    parts = eid.split(':', 1)
    if parts[0] in ["gene","transcript","CDS"]:
	return parts[1]
    return eid

def feature(f):
    attrs = f[8]
    if 'ID' in attrs:
        attrs['ID'] = stripPrefix(attrs['ID'])
    if 'Parent' in attrs:
        attrs['Parent'] = map(stripPrefix, attrs['Parent'])
    return f
