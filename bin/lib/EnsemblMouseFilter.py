import re
from .GffFilter import GffFilter
from urllib.request import urlopen

class EnsemblMouseFilter (GffFilter) :
    # index mapping ensembl IDs to MGI ids
    # Shared by all instances. The first one to try to access the index creates it.
    EID2MGI = None
    def getEid2MgiIndex (self) :
        if self.EID2MGI:
            return self.EID2MGI
        self.EID2MGI = {}
        # URL for a query that returns all mouse genes with their Ensembl IDs. Data returned 
        # as TSV with 3 columns: MGIid, symbol, Ensembl id
        url = '''https://www.mousemine.org/mousemine/service/query/results?query=%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22Gene.primaryIdentifier+Gene.symbol+Gene.crossReferences.identifier%22+longDescription%3D%22%22+sortOrder%3D%22Gene.crossReferences.identifier+asc%22+constraintLogic%3D%22A+and+B%22%3E%3Cconstraint+path%3D%22Gene.crossReferences.source.name%22+code%3D%22A%22+op%3D%22%3D%22+value%3D%22Ensembl+Gene+Model%22%2F%3E%3Cconstraint+path%3D%22Gene.dataSets.name%22+code%3D%22B%22+op%3D%22%3D%22+value%3D%22Mouse+Gene+Catalog+from+MGI%22%2F%3E%3C%2Fquery%3E&format=tab''' 

        self.log("Getting MGI/Ensembl ID associations from: " + url)
        for r in urlopen(url):
            rr = r.decode('utf-8').strip().split()
            self.EID2MGI[rr[2]] = rr
        return self.EID2MGI

    def stripPrefix(self, eid) :
        parts = eid.split(':', 1)
        if parts[0] in ["gene","transcript","CDS"]:
            return parts[1]
        return eid

    def processFeature(self, feat):
        attrs = feat[8]
        if feat[2] == "gene" and attrs.get('biotype', None) == 'protein_coding':
            feat[2] = 'protein_coding_gene'
        if 'ID' in attrs:
            attrs['ID'] = self.stripPrefix(attrs['ID'])
        if 'Parent' in attrs:
            attrs['Parent'] = list(map(self.stripPrefix, attrs['Parent']))
        if "projection_parent_gene" in attrs:
            eid = attrs['projection_parent_gene'].split(".")[0]
            mgi = self.getEid2MgiIndex().get(eid, None)
            if mgi:
                attrs['curie'] = mgi[0]
                attrs['Name'] = mgi[1]
            else:
                attrs['curie'] = eid
                attrs['Name'] = eid
        if "transcript_id" in attrs:
            attrs["transcript_id"] = "ENSEMBL:" + attrs["transcript_id"]
        if "protein_id" in attrs:
            attrs["protein_id"] = "ENSEMBL:" + attrs["protein_id"]
        return feat

