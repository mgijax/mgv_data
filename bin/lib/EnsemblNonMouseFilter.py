from .GffFilter import GffFilter

class EnsemblNonMouseFilter (GffFilter) :
    # used for finding/extracting canonical id from another attribute, by taxon id
    taxon2re = {
        "7955":  ('description', re.compile(r'Source:ZFIN.*Acc:([A-Z0-9-]+)'), "ZFIN:"),
        "9606":  ('description', re.compile(r'Acc:HGNC:(\d+)'),                "HGNC:"),
        "10116": ('description', re.compile(r'Source:RGD.*Acc:(\d+)'),         "RGD:"),
        "559292":('description', re.compile(r'Source:SGD.*Acc:(S\d+)'),        "SGD:"),
        "7227":  ('gene_id',     re.compile(r'(.*)'),                          "FB:"), 
        "6239":  ('gene_id',     re.compile(r'(.*)'),                          "WB:"),

    }
    def __init__(self, *args):
        GffFilter.__init__(self, *args)
        self.attrName, self.id_re, self.id_prefix = self.taxon2re[GCONFIG["taxonid"]]

    def processFeature (self, f) :
        attrs = f[8]
        if f[2] == "gene" and attrs.get('biotype', None) == 'protein_coding':
            f[2] = 'protein_coding_gene'
        attrVal = attrs.get(self.attrName,'')
        m = self.id_re.search(attrVal)
        if m:
            attrs['curie'] = self.id_prefix + m.group(1)
        return f
        
