import re
from .Filter import Filter
from .gff3lite import parseLine, formatLine, Gff3Parser

class GffFilter (Filter) :
    def __init__(self, src, gcfg, dcfg):
        src = Gff3Parser(src, returnHeader=True, returnGroups=True).iterate()
        Filter.__init__(self, src, gcfg, dcfg)
        self.currChr = None

    def processNext (self, obj) :
        if type(obj) is str:
            # the header
            return obj
        else:
            # obj is a list of features
            m = self.processModel(obj)
            return ''.join(m) if m else None

    def processModel(self, model):
        model = list(filter(lambda x: x, [self._processFeature(f) for f in model]))
        return model

    def _processFeature (self, f):
        # universal transform: filter out features on non-matching chromosome 
        if not self.chr_re.match(f[0]):
            return None
        # universal transform: filter out features with excluded types
        if 'include_types' in self.dcfg and f[2] not in self.dcfg['include_types']:
            return None
        if 'exclude_types' in self.dcfg and f[2] in self.dcfg['exclude_types']:
            return None
        if f[0] != self.currChr:
            self.log("Chromosome " + f[0])
            self.currChr = f[0]
        ff = self.processFeature(f)
        return formatLine(ff) if ff else None

CURIE_INFO = [{
  "prefix" : "ENSEMBL",
  "baseRegex" : re.compile(r'ENS[A-Z]+[GTP]+\d+'),
}, {
  "prefix" : "ENSEMBL",
  "baseRegex" : re.compile(r'MGP_[A-Z0-9]+_[GTP]\d+'),
}, {
  "prefix" : "RefSeq",
  "baseRegex" : re.compile(r'[NX][RMP]_\d+')
}]

def curie_ize (ident) :
    i = ident.find(":")
    if i >= 0:
        return ident
    for info in CURIE_INFO:
        if info["baseRegex"].match(ident):
            return "%s:%s" % (info["prefix"], ident)
    return None

