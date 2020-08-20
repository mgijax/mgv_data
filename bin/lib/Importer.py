import re
import os
import os.path
import json
import sys
import gzip
from .Filter import filterNameMap
from .gff3lite import Gff3Parser, parseLine, formatLine


class Importer :
    def __init__ (self, builder, type, cfg, output_dir, debug) :
        self.builder = builder
        self.debug = debug
        self.type = type
        self.cfg = cfg
        tcfg = self.cfg[self.type]
        self.chr_re = self.cfg["chr_re"]
        self.exclude_types = tcfg.get("exclude_types", [])
        self.fpath = tcfg["fpath"]
        self.log = self.builder.log
        self.output_dir = output_dir
        self.filters = tcfg.get("filters", [])[:]

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

    # OVERRIDE ME
    def processLine(self, line):
        pass

    def go (self) :
        self.log("Importing file: " + self.fpath)
        self.log("Filters: " + str(self.filters))
        count = 0
        if self.debug:
            return
        for line in self.filterDownloadedFile():
            self.processLine(line)

class FastaImporter (Importer) :
    def go (self) :
        lfunc = lambda s: self.log(s, newline='')
        odir = os.path.join(self.output_dir, self.cfg["name"], "assembly")
        self.log("Importing Fasta: %s -> %s/assembly" % (self.fpath, odir))
        if not self.debug:
            self.builder.ensureDirectory(odir, empty=True)
            self.split(self.filterDownloadedFile(), odir, self.chr_re, lfunc)

    def split (self, ifd, odir, chr_re, logFcn):
      ofd = None
      writing = False
      lcount = 0
      try:
          line = next(ifd)
          while line:
              lcount += 1
              if line.startswith('>'):
                  seqid = line.split()[0][1:]
                  writing = chr_re.match(seqid)
                  if writing:
                      # output is plain text, not Fasta
                      ofile = os.path.join(odir, seqid + '.txt')
                      if ofd: ofd.close()
                      ofd = open(ofile, 'w')
                      logFcn(line)
                  else:
                      logFcn("Skipping: " + line)
                  line = next(ifd)
              else:
                  if writing:
                      ofd.write(line[:-1])
                  line = next(ifd)
      except StopIteration:
          pass
      if ofd: ofd.close()

class GffImporter (Importer) :
    def __init__ (self, *args):
        Importer.__init__(self, *args)
        self.doSort = self.cfg[self.type].get("doSort", False)
        self.transcriptChunkSize = int(self.cfg[self.type]["transcriptChunkSize"])

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
        return gffStream

    def filterDownloadedFile (self) :
        for obj in Importer.filterDownloadedFile(self):
            if type(obj) is str:
                # yield the header
                yield obj
            else:
                # yield features on matching chromosomes
                if len(obj):
                    if self.chr_re.match(obj[0][0]):
                        yield obj

    def getFileName(self, track, chr, blk):
        if chr:
            return os.path.join(self.mOutputDir, track, chr, str(blk) + '.gff3')
        else:
            return os.path.join(self.mOutputDir, track, str(blk), '.gff3')

    def getFileHandle(self, fname):
        if fname == self.currFileName:
            # currently open file
            return self.currFile
        elif fname in self.outputFiles:
            # different from current but seen before, reopen in append mode
            self.currFile.close()
            self.currFileName = fname
            self.currFile = open(fname, 'a')
            return self.currFile
        else:
            # first time seeing this fname
            self.builder.ensureDirectory(os.path.dirname(fname))
            self.currFile and self.currFile.close()
            self.currFileName = fname
            self.currFile = open(fname, 'w')
            self.outputFiles[fname] = 0

            fn = os.path.basename(fname)
            dn = os.path.basename(os.path.dirname(fname))
            self.log('%s/%s ' % (dn, fn), newline=' ', timestamp=False)
            return self.currFile

    def formatFeature(self, f):
        return formatLine(f)

    def writeGrpToBlk(self, grp, track, chr, blk):
        fname = self.getFileName(track, chr, blk)
        fd = self.getFileHandle(fname)
        for f in grp:
            if 'exons' in f[8]:
                fd.write(self.formatFeature(f))
                self.outputFiles[fname] += 1

    def writeGrp (self, track, grp) :
        if len(grp) == 0: return
        chr = grp[0][0]
        if self.transcriptChunkSize == 0 :
            self.writeGrpToBlk(grp, track, None, 0)
        elif self.transcriptChunkSize == 1:
            self.writeGrpToBlk(grp, track, chr, 0)
        else:
            gStart = min([f[3] for f in grp])
            gEnd = max([f[4] for f in grp])
            startBlk = gStart // self.transcriptChunkSize
            endBlk = gEnd // self.transcriptChunkSize
            for blk in range(startBlk, endBlk+1):
                self.writeGrpToBlk(grp, track, chr, blk)

    def writeGene(self, f):
        self.tlFile.write(self.formatFeature(f))

    # Processes a gene model group.
    # For transcripts, adds an exons attribute whose value is a list 
    # of exon coordinates. Each exon is encoded as 'offset_length', where offset is equal to
    # exon.start - transcript.start, and length is equal to exon.end - exon.start + 1.
    # Also keeps track of all chromosomes seen in the file and the max coordinate seen on each.
    # Writes genes to genes track and transcripts to transcripts track.
    def processGrp(self, grp):
        # Build index from ID to the features that have a Parent attribute with that ID.
        # Keep track of the chromosome in the file, their order, and the max coordinate of features on each.
        pid2kids = {}
        for f in grp:
          seqid = f[0]
          c = self.seenChrs.get(seqid, None)
          if not c:
            c = { 'name': seqid, 'length': 0 }
            self.seenChrs[seqid] = c
            self.seenChrOrder.append(c)
          c['length'] = max(c['length'], f[4])
          for pid in f[8].get('Parent',[]):
            pid2kids.setdefault(pid, []).append(f)
        # For transcripts, encode exon coordinates and store in exons attribute
        transcripts = []
        for f in grp:
          fid = f[8].get('ID', None)
          if 'Parent' not in f[8]:
            # top level feature
            self.writeGene(f)
          elif fid in pid2kids:
            # mid level feature
            # For efficiency of transfer, encode the exons for a transcript into a col9 attribute.
            # For each exon, just need its offset from the start of the transcript and its length.
            # This requires a lot less space (~time) to transfor than representing each exon as a separate
            # feature with full coordinates. NOTE that we lose exon IDs (and potentially other attributes) 
            # by doing this. Maybe worth it, maybe not?
            exons = list([x for x in pid2kids[fid] if x[2] == 'exon'])
            exons.sort(key=lambda e: e[3])
            eexons = f[8].setdefault('exons',[])
            for e in exons:
              offset = e[3] - f[3]
              length = e[4] - e[3] + 1
              eexons.append('%d_%d' % (offset,length))
            #
            cdss = list([x for x in pid2kids[fid] if x[2] == 'CDS'])
            if cdss:
              cdss.sort(key=lambda x: x[3])
              attrs = cdss[0][8]
              cid = attrs['ID']
              prid = attrs.get('protein_id', attrs.get('curie', ''))
              cstart = cdss[0][3]
              cend = cdss[-1][4]
              f[8]['cds'] = '%s|%s|%d|%d' % (cid, prid, cstart, cend)
            #
            transcripts.append(f)
          else:
            # leaf level feature
            pass
        self.writeGrp('transcripts', transcripts)

    def go (self):
        self.mOutputDir = os.path.join(self.output_dir, self.cfg["name"], "models")
        self.gOutputDir = os.path.join(self.mOutputDir, "genes")
        self.tOutputDir = os.path.join(self.mOutputDir, "transcripts")
        self.builder.ensureDirectory(self.gOutputDir, empty=True)
        self.builder.ensureDirectory(self.tOutputDir, empty=True)
        #
        self.seenChrs = {}
        self.seenChrOrder = []
        self.tlFileName = os.path.join(self.gOutputDir, '0.gff3')
        self.tlFile = open(self.tlFileName, 'w')
        self.log("Writing genes to: " + self.tlFileName)
        self.currFileName = None
        self.currFile = None
        self.outputFiles = {}
        #
        self.log("Importing GFF: %s -> %s" % (self.fpath, self.mOutputDir))
        if not self.debug:
            for grp in self.filterDownloadedFile():
                if type(grp) is not str:
                    self.processGrp(grp)
##
class OrthologyImporter (Importer) :

    def __init__(self, *args):
        Importer.__init__(self, *args)
        self.taxon2file = {}
        self.inCount = 0
        self.output_dir = os.path.join(self.output_dir, self.cfg["name"], "orthology")
        self.builder.ensureDirectory(self.output_dir)

    #
    def getFile (self, taxonid) :
      fname = os.path.join(self.output_dir, taxonid + '.json')
      if taxonid not in self.taxon2file:
        self.taxon2file[taxonid] = {
          "fname" : fname,
          "fd" : open(fname, 'w'),
          "count" : 0
        }
      return self.taxon2file[taxonid]
      
    #
    def closeAll (self) :
      for rec in list(self.taxon2file.values()):
        rec["fd"].write("]")
        rec["fd"].close()
    #
    def outputRecord (self, id1, tx1, id2, tx2, yn) :
      rec = self.getFile(tx1)
      if rec["count"] == 0:
          rec["fd"].write("[")
      else:
          rec["fd"].write(",")
      rec["count"] += 1
      rec["fd"].write(json.dumps([id1, tx1, id2, tx2, yn]))
      rec["fd"].write("\n")

    def processLine (self, line) :
      if line.startswith("#"):
        return
      # skip column labels
      self.inCount += 1
      if self.inCount == 1:
        return
      # parse and ouput
      fs = line[:-1].split('\t')
      # Example desired record - just the essentials:
      # ['FB:FBgn0033981', 'NCBITaxon:7227', 'MGI:3646373', 'NCBITaxon:10090', 'YN']
      self.outputRecord(
        fs[0],
        fs[2].replace('NCBITaxon:',''),
        fs[4],
        fs[6].replace('NCBITaxon:',''),
        fs[11][0] + fs[12][0])

    def go (self) :
        Importer.go(self)
        self.closeAll()

importerNameMap = {
    "models" : GffImporter,
    "assembly" : FastaImporter,
    "orthologs" : OrthologyImporter,
}
