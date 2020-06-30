import re
#[Source:RGD Symbol;Acc:1311118]
RE = re.compile(r'Source:RGD.*Acc:(\d+)')
def feature(f):
    attrs = f[8]
    descr = attrs.get('description','')
    m = RE.search(descr)
    if m:
        attrs['cID'] = 'RGD:'+m.group(1)
    return f
