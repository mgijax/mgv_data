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
import json


class Filter:
    def __init__ (self, src, gcfg, dcfg):
        self.src = src
        self.gcfg = gcfg
        self.dcfg = dcfg
        self.chr_re = re.compile(self.dcfg.get("chr_re", ".*"))

    def log (self, s) :
        sys.stderr.write(s)
        sys.stderr.write('\n')
        sys.stderr.flush()

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

