import re
RE = re.compile(r'Acc:(HGNC:\d+)')
def feature(f):
    attrs = f[8]
    descr = attrs.get('description','')
    m = RE.search(descr)
    if m:
        attrs['cID'] = m.group(1)
    return f
