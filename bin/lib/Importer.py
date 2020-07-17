import re
import os
import sys
import gzip
from .Filter import filterNameMap
from .gff3lite import Gff3Parser, parseLine, formatLine

class Importer :
    def __init__ (self, builder, type, cfg, output_dir) :
        self.builder = builder
        self.type = type
        self.cfg = cfg
        tcfg = self.cfg[self.type]
        self.chr_re = self.cfg["chr_re"]
        self.exclude_types = tcfg.get("exclude_types", [])
        self.fpath = tcfg["fpath"]
        self.log = self.builder.log
        self.output_dir = output_dir
        self.filters = tcfg.get("filters", [])[:]

    def ensureDirectory (self, d):
        if not os.path.exists(d):
            os.makedirs(d)

    def streamDownloadedFile (self) :
        if self.fpath.endswith(".gz") :
            fp = gzip.open(self.fpath)
            decode = True
        else:
            fp = open(self.fpath)
            decode = False
        #
        for line in fp:
            if decode:
                line = line.decode('utf-8')
            yield line
        fp.close()
        
    def filterDownloadedFile (self) :
        src = self.streamDownloadedFile()
        if hasattr(self, "filters"):
            for filt in self.filters:
                filtCls = filterNameMap[filt]
                src = filtCls(self)(src)
        return src

    def go (self) :
        self.log("Importing file: " + self.fpath)
        self.log("Filters: " + str(self.filters))
        count = 0
        for line in self.filterDownloadedFile():
            self.log(line, timestamp=False, newline='')
            count += 1
            if count == 10:
                break

from .FastaImporter import split
class FastaImporter (Importer) :
    def go (self) :
        odir = os.path.join(self.output_dir, self.cfg["name"], "sequences")
        self.ensureDirectory(odir)
        split(self.filterDownloadedFile(), odir, self.chr_re, self.log)

from .GffImporter import Gff3Importer, getArgs
class GffImporter (Importer) :
    def __init__ (self, *args):
        Importer.__init__(self, *args)
        self.doSort = self.cfg[self.type].get("doSort", False)

    def streamDownloadedFile (self) :
        # GffImporter streams the file a gene model at a time.
        # Each yielded item is a list of features, except for the first item, which is the header.
        #
        # First, the thing that streams the file one line at a time
        src = Importer.streamDownloadedFile(self)
        # Pass to a Gff3Parser, set to return models
        gp = Gff3Parser(src, returnHeader=True, returnGroups=True)
        # File must be sorted properly. doSort indicates that a sort is needed.
        # This is an INTERNAL sort, ie, this loads the whole file into memory.
        # FIXME: would be more efficient to sort it externally, and store it that way.
        if self.doSort:
            gffStream = gp.sortIterate()
        else:
            gffStream = gp.iterate()
        for x in gffStream:
            yield x

    def filterDownloadedFile (self) :
        for obj in Importer.filterDownloadedFile(self):
            if type(obj) is str:
                yield obj
            else:
                for f in obj:
                    yield formatLine(f)

    def go (self):
        opts = getArgs([
            "-p", self.cfg["name"],
            "-g", self.cfg["label"],
            "-x", self.cfg["taxonid"],
            "-d", self.output_dir,
            "-k", "4000000"
        ])
        imp = Gff3Importer(self.filterDownloadedFile(), opts)
        imp.main()
