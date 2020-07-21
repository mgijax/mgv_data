
import sys
import re
from .gff3lite import parseLine, formatLine
from urllib.request import urlopen

class Filter:
    def __init__ (self, impobj):
        self.importer = impobj
        self.log = self.importer.log

    # OVERRIDE me
    def processObj(self, obj):
        return obj

    def __call__(self, src):
        for obj in src:
            obj = self.processObj(obj)
            if obj:
                yield obj

class GffFilter (Filter) :
    def processObj (self, obj) :
        if type(obj) is str:
            # the header
            return obj
        else:
            # obj is a list of features
            return self.processModel(obj)

    def processModel(self, model):
        return list(filter(lambda x: x, [self.processFeature(f) for f in model]))


class EnsemblMouseStrain (GffFilter) :
    # index mapping ensembl IDs to MGI ids
    # Shared by all instances. The first one to try to access the index creates it.
    EID2MGI = None
    def getEid2MgiIndex (self) :
        if self.EID2MGI:
            return self.EID2MGI
        self.EID2MGI = {}
        # URL for a query that returns all mouse genes with their Ensembl IDs. Data returned 
        # as TSV with 3 columns: MGIid, symbol, Ensembl id
        url = '''http://www.mousemine.org/mousemine/service/query/results?query=%3Cquery+name%3D%22%22+model%3D%22genomic%22+view%3D%22Gene.primaryIdentifier+Gene.symbol+Gene.crossReferences.identifier%22+longDescription%3D%22%22+sortOrder%3D%22Gene.crossReferences.identifier+asc%22+constraintLogic%3D%22A+and+B%22%3E%3Cconstraint+path%3D%22Gene.crossReferences.source.name%22+code%3D%22A%22+op%3D%22%3D%22+value%3D%22Ensembl+Gene+Model%22%2F%3E%3Cconstraint+path%3D%22Gene.dataSets.name%22+code%3D%22B%22+op%3D%22%3D%22+value%3D%22Mouse+Gene+Catalog+from+MGI%22%2F%3E%3C%2Fquery%3E&format=tab'''

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
        if not self.importer.chr_re.match(feat[0]):
            return None
        if feat[2] in self.importer.exclude_types:
            return None
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
                attrs['cID'] = mgi[0]
                attrs['Name'] = mgi[1]
        return feat

class EnsemblNonMouseFilter (GffFilter) :
    taxon2re = {
        "7955":  ('description', re.compile(r'Source:ZFIN.*Acc:([A-Z0-9-]+)'), "ZFIN:"),
        "9606":  ('description', re.compile(r'Acc:(HGNC:\d+)'),                ""),
        "10116": ('description', re.compile(r'Source:RGD.*Acc:(\d+)'),         "RGD:"),
        "559292":('description', re.compile(r'Source:SGD.*Acc:(S\d+)'),        "SGD:"),
        "7227":  ('gene_id',     re.compile(r'(.*)'),                          "FB:"), 
        "6239":  ('gene_id',     re.compile(r'(.*)'),                          "WB:"),

    }
    def __init__(self, *args):
        GffFilter.__init__(self, *args)
        self.cfg = self.importer.cfg
        self.attrName, self.id_re, self.id_prefix = self.taxon2re[self.cfg["taxonid"]]

    def processFeature (self, f) :
        attrs = f[8]
        attrVal = attrs.get(self.attrName,'')
        m = self.id_re.search(attrVal)
        if m:
            attrs['cID'] = self.id_prefix + m.group(1)
        return f
        
class MgiGff (GffFilter) :
    def __init__ (self, impobj) :
        Filter.__init__(self, impobj)
        self.i2i = {}

    def processFeature(self, f):
        attrs = f[8]
        if f[2] in ["gene","pseudogene"]:
            #
            f[2] = attrs['so_term_name']
            del attrs['so_term_name']
            #
            attrs['cID'] = attrs['curie']
            del attrs['curie']
        elif f[2] == "CDS":
            pid = attrs.get('protein_id', attrs['ID'])
            attrs['ID'] = pid
        elif f[2] != "exon":
            tid = attrs.get('transcript_id', attrs['ID'])
            self.i2i[attrs['ID']] = tid
            attrs['ID'] = tid
        if 'Parent' in attrs:
            attrs['Parent'] = [self.i2i.get(p, p) for p in attrs['Parent']]
        return f

class AllianceGff (GffFilter) : 
    def processFeature(self, f) :
        attrs = f[8]
        attrs.pop("description", None)
        if "curie" in attrs and "Parent" not in attrs:
            attrs["cID"] = attrs["curie"]
        return f

class SgdGff (AllianceGff) : 
    # SGF gff issues:
    # Yeast transcription/transpation is simpler b.c. no introns. Just an mRNA and a CDS.
    # Therefore in the GFF:
    # - CDS features do not have IDs (no need to tie multiple CDSs together)
    #   => assign IDs based on Parent's ID
    # - CDS features can preceed their mRNA parents (forward reference issue)
    #   => have to re-sort the model's features
    def processModel(self, model) :
        # re-sort by: level, then start pos. 
        def keyfun(f):
            l1 = 0
            if "Parent" in f[8]:
                l1 = 1 if f[2] in ["transcript", "mRNA"] else 2
            return (l1, f[4])
        model.sort(key=keyfun)
        exons = []
        for f in model:
            if f[2] in ["transcript","mRNA"]:
               e = f[:]
               e[2] = "exon"
               e[8] = { "Parent" : f[8]["ID"] }
               exons.append(e)
            if f[2] == "CDS":
                pid = f[8]["Parent"][0]
                cid = pid.replace("mRNA", "CDS")
                if cid == pid:
                    cid = pid + "_CDS"
                f[8]["ID"] = [cid]
        model = model + exons
        return AllianceGff.processModel(self, model)

#
filterNameMap = {
  "ensemblMouseStrain" : EnsemblMouseStrain,
  "ensemblNonMouseFilter" : EnsemblNonMouseFilter,
  "mgiGff" : MgiGff,
  "allianceGff" : AllianceGff,
  "sgdGff" : SgdGff
}
