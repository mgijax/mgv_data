
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

# URL for a query that returns all mouse genes with their Ensembl IDs. Data returned 
# as TSV with 3 columns: MGIid, symbol, Ensembl id
url = '''http://www.mousemine.org/mousemine/service/query/results?query=%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22Gene.primaryIdentifier+Gene.symbol+Gene.crossReferences.identifier%22+longDescription%3D%22%22+sortOrder%3D%22Gene.crossReferences.identifier+asc%22+constraintLogic%3D%22A+and+B%22%3E%3Cconstraint+path%3D%22Gene.crossReferences.source.name%22+code%3D%22A%22+op%3D%22%3D%22+value%3D%22Ensembl+Gene+Model%22%2F%3E%3Cconstraint+path%3D%22Gene.dataSets.name%22+code%3D%22B%22+op%3D%22%3D%22+value%3D%22Mouse+Gene+Catalog+from+MGI%22%2F%3E%3C%2Fquery%3E&format=tab'''

eid2mgi = {}
for r in urlopen(url):
    rr = r.decode('utf-8').strip().split()
    eid2mgi[rr[2]] = rr

def feature(f):
    attrs = f[8]
    if "projection_parent_gene" in attrs:
        eid = attrs['projection_parent_gene'].split(".")[0]
        mgi = eid2mgi.get(eid, None)
        if mgi:
            attrs['cID'] = mgi[0]
            attrs['Name'] = mgi[1]
    return f
