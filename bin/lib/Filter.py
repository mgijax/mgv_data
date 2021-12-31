#
# Filter.py
#
# Herein lies all the ugliness needed to convert GFF3 files as downloaded from their
# respective sources into GFF3 as needed for import. 
# A filter implements a set of transformations on a gene model.
# In genomes.json, each stanza specifies (by name) zero or more filters to apply during import.
# The names used in genome.json are defined in filterNameMap (at the end of this file), which maps
# each name to the class that implements it.
# A filter class is a subclass of GffFilter that defines one or both methods
# processModel and processFeature. 
# Each model is a list of individual features (gene, transcripts, exons, CDSs, etc)
# that make up the model. The general flow is that processModel is called with the
# list of features, and it then calls processFeature on each list element in turn.
# The top level feature (the gene or pseudogene) is always the first
# element of the list. 
#
# Overriding processFeature:
# - good for simple line-by-line transforms
# 

import sys
import os
import re
from .gff3lite import parseLine, formatLine, Gff3Parser
import json
from urllib.request import urlopen

GCONFIG = json.loads(os.environ.get("GCONFIG", "{}"))
DCONFIG = json.loads(os.environ.get("DCONFIG", "{}"))
sys.stderr.write('GCONFIG=' + str(GCONFIG) + '\n')
sys.stderr.write('DCONFIG=' + str(DCONFIG) + '\n')

CHR_RE = re.compile(DCONFIG.get("chr_re", ".*"))

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

class Filter:
    def __init__ (self, src):
        self.src = src

    def log (self, s) :
        sys.stderr.write(s)
        sys.stderr.write('\n')

    def __iter__(self):
        return self

    def __next__ (self):
        res = None
        while res is None:
            obj = next(self.src)
            res = self.processNext(obj)
        return res

    def processNext(self, obj):
        return obj

