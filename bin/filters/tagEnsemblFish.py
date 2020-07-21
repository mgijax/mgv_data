import re
#[Source:ZFIN%3BAcc:ZDB-GENE-020419-25]
RE = re.compile(r'Source:ZFIN.*Acc:([A-Z0-9-]+)')
def feature(f):
    attrs = f[8]
    descr = attrs.get('description','')
    m = RE.search(descr)
    if m:
        attrs['cID'] = 'ZFIN:'+m.group(1)
    return f

