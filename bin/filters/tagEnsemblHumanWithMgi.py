import re
try:
    from urllib.request import urlopen
except:
    from urllib import urlopen

# URL for a query that returns all human genes and their mouse orthologs.
# 2 columns: HGNC id, MGI id
url="http://www.mousemine.org/mousemine/service/query/results?query=%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22Gene.crossReferences.identifier+Gene.homologues.homologue.primaryIdentifier%22+longDescription%3D%22%22+constraintLogic%3D%22A+and+B+and+C%22%3E%3Cconstraint+path%3D%22Gene.crossReferences.source.name%22+code%3D%22A%22+op%3D%22%3D%22+value%3D%22HGNC%22%2F%3E%3Cconstraint+path%3D%22Gene.organism.taxonId%22+code%3D%22B%22+op%3D%22%3D%22+value%3D%229606%22%2F%3E%3Cconstraint+path%3D%22Gene.homologues.homologue.organism.taxonId%22+code%3D%22C%22+op%3D%22%3D%22+value%3D%2210090%22%2F%3E%3C%2Fquery%3E&format=tab"

hid2mgi = {}
for r in urlopen(url):
    rr = r.decode('utf-8').strip().split()
    hid2mgi[rr[0]] = rr[1]

RE = re.compile(r'Acc:(HGNC:\d+)')
def feature(f):
    attrs = f[8]
    descr = attrs.get('description','')
    m = RE.search(descr)
    if m:
 	mgiid = hid2mgi.get(m.group(1), '')
	if mgiid:
	    attrs['cID'] = mgiid
    return f
