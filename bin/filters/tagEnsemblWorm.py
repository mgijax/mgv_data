#gene_id=WBGene00022278
def feature(f):
    attrs = f[8]
    gid = attrs.get('gene_id','')
    if gid:
        attrs['cID'] = 'WB:'+gid
    return f
