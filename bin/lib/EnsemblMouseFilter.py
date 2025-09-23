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

        # URL for an MGI database report of gene models. The columns we need are 0 (MGI id), 1 (symbol)
        # and 10 (ENSEMBL gene id).
        url = '''https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt'''

        self.log("Getting MGI/Ensembl ID associations from: " + url)
        first = True
        for r in urlopen(url):
            # skip the header line
            if first:
                first = False
                continue
            rr = r.decode('utf-8').strip().split('\t')
            if rr[10] != 'null':
                self.EID2MGI[rr[10]] = [rr[0], rr[2], rr[10]]
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

