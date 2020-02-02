#gene_id=FBgn0002121
def feature(f):
    attrs = f[8]
    gid = attrs.get('gene_id','')
    if gid:
 	attrs['cID'] = 'FB:'+gid
    return f
