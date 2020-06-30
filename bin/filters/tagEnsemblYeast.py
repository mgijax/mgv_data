import re
#[Source:SGD%3BAcc:S000028594]
RE = re.compile(r'Source:SGD.*Acc:(S\d+)')
def feature(f):
    attrs = f[8]
    descr = attrs.get('description','')
    m = RE.search(descr)
    if m:
        attrs['cID'] = 'SGD:'+m.group(1)
    return f
