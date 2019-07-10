#
#
# importGff3.py
#
# usage:
#    python importGff3.py < MyGenomeAnnotations.gff3
#
# Chunking:
#    For efficiency, output can be chunked at specified granularities.
#       chunkSize   meaning
#       0           No chunking. Everything in one file. The file's name is the track name.
#                   It is located in the root directory for the genome.
#       1           Chunk by chromosome. Track name = subdirectory in genome root dir.
#                   Then one subdir per chromosome below this.One file per chromosome.
#       > 1         Chunk size. As above, but output per chromosome is chunked by this size block factor.
# 
import types
import os
import sys
from gff3lite import Gff3Parser, formatLine
import argparse
import json
import time

TIMESTAMP = time.asctime(time.localtime(time.time()))

#
class Gff3Importer:
  def __init__(self, gffSource, opts):
    self.opts = opts
    self.parser = Gff3Parser(gffSource, returnHeader=True, returnGroups=True, convertDots='.')
    self.datastream = self.parser.iterate()
    self.genomeInfo = {
      "name" : None,
      "timestamp" : None,
      "chromosomes" : None,
      "tracks" : None
    }
    #
    self.seenChrs = {}
    self.seenChrOrder = []
    #
    self.outputDir = None
    self.outputFiles = {}
    self.currFileName = None
    self.currFile = None
    # top-level file name and file handle
    self.tlFileName = None
    self.tlFile = None
    #
    self.processHeader()

  def getOpt(self, n):
    v = getattr(self.opts, n, None)
    if v is None:
      return None
    #
    if v.startswith('##'):
      v = self.pragmas[v[2:]]
    return v

  def processHeader (self) :

    self.pragmas = next(self.datastream)
    self.genomeInfo['taxonid'] = self.getOpt('taxonid')
    self.genomeInfo['name'] = self.getOpt('genome')
    self.genomeInfo['timestamp'] = self.getOpt('timestamp')
    #
    chrs = []
    copt = self.getOpt('chromosomes')
    if copt:
      if type(copt) is str:
        # parse the command line arg
        lst = map(lambda x: x.split(':'), copt.split(','))
        for item in lst:
          chrs.append({
            "name": item[0],
            "length": int(item[-1]) if len(item) > 1 else 0
          })
      else:
        # parse the 'sequence-region' pragma lines.
        for c in copt:
          item = c.split()
          chrs.append({
            "name": item[0],
            "length": int(item[-1])
          })
      self.genomeInfo['chromosomes'] = chrs
    #
    self.genomeInfo['tracks'] = [{
      'name': 'genes',
      'reader' : {
        'type' : 'ChunkedGff3FileReader',
        'chunkSize' : 0
      }
    }, {
      'name': 'transcripts',
      'reader' : {
        'type' : 'ChunkedGff3FileReader',
        'chunkSize': self.opts.chunkSize
      }
    }, {
      'name': 'sequences',
      'reader' : {
        'type' : 'MouseMineSequenceReader',
        'url' : 'http://www.mousemine.org/mousemine/service'
      }
    }]
    self.outputDir = os.path.join(self.opts.outputDir, self.sanitizeName(self.genomeInfo['name']))
    
  def log(self, s):
    sys.stderr.write(s)

  def sanitizeName (self, n):
    return n.replace('/','').lower()

  def ensureDirectory (self, d):
    if not os.path.exists(d):
      os.makedirs(d)

  def getFileName(self, track, chr, blk):
    if chr:
      return os.path.join(self.outputDir, track, chr, str(blk))
    else:
      return os.path.join(self.outputDir, track, str(blk))

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
      self.ensureDirectory(os.path.dirname(fname))
      self.currFile and self.currFile.close()
      self.currFileName = fname
      self.currFile = open(fname, 'w')
      self.outputFiles[fname] = 0
      self.log(fname + '\n')
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
    if self.opts.chunkSize == 0 :
      self.writeGrpToBlk(grp, track, None, 0)
    elif self.opts.chunkSize == 1:
      self.writeGrpToBlk(grp, track, chr, 0)
    else:
      gStart = min(map(lambda f: f[3], grp))
      gEnd = max(map(lambda f: f[4], grp))
      startBlk = gStart // self.opts.chunkSize
      endBlk = gEnd // self.opts.chunkSize
      for blk in range(startBlk, endBlk+1):
        self.writeGrpToBlk(grp, track, chr, blk)

  def writeGene(self, f):
    self.tlFile.write(self.formatFeature(f))

  # Processes a gene model group. For top level features, adds a tCount attribute that contains the
  # number of transcripts for that gene. For transcripts, adds an exons attribute whose value is a list 
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
    # For genes, count transcripts and store in tCount attribute.
    # For transcripts, encode exon coordinates and store in exons attribute
    transcripts = []
    for f in grp:
      fid = f[8].get('ID', None)
      if 'Parent' not in f[8]:
        # top level feature
        f[8]['tCount'] = len(pid2kids.get(fid,[]))
        self.writeGene(f)
      elif fid in pid2kids:
        # mid level feature
        exons = list(filter(lambda x: x[2] == 'exon', pid2kids[fid]))
        exons.sort(key=lambda e: e[3])
        eexons = f[8].setdefault('exons',[])
        for e in exons:
          offset = e[3] - f[3]
          length = e[4] - e[3] + 1
          eexons.append('%d_%d' % (offset,length))
        #
        cdss = list(filter(lambda x: x[2] == 'CDS', pid2kids[fid]))
        if cdss:
          cdss.sort(key=lambda x: x[3])
          cid = cdss[0][8]['ID']
          cstart = cdss[0][3]
          cend = cdss[-1][4]
          f[8]['cds'] = '%s|%d|%d' % (cid, cstart, cend)
        #
        transcripts.append(f)
      else:
        # leaf level feature
        pass
    self.writeGrp('transcripts', transcripts)

  def writeGenomeInfo(self):
    fd = open(os.path.join(self.outputDir, 'index.json'), 'w')
    fd.write(json.dumps(self.genomeInfo, indent=2))
    fd.close()

  def main (self) :
    self.ensureDirectory(self.outputDir)
    self.ensureDirectory(os.path.join(self.outputDir, 'genes'))
    self.tlFileName = os.path.join(self.outputDir, 'genes', '0')
    self.tlFile = open(self.tlFileName, 'w')
    for i,grp in enumerate(self.datastream):
      toplevel = self.processGrp(grp)
      if self.opts.sample and i > 1000:
         break
    #
    if self.genomeInfo['chromosomes'] is None:
      self.genomeInfo['chromosomes'] = self.seenChrOrder
    #
    self.writeGenomeInfo()
#
def getArgs () :
  parser = argparse.ArgumentParser(description='Import genome annotations from a GFF3 file.')
  parser.add_argument('-g','--genomeName',
    metavar='NAME',
    default='##genome-name',
    dest='genome',
    help='Name of the genome. Required. Either specify the name directly (eg, "C57BL/6J") or provide the name of a pragma in the GFF3 header (eg "##genome-name")')
  parser.add_argument('-x','--taxonid',
    metavar='TAXONID',
    default='##taxonid',
    dest='taxonid',
    help='NCBI taxon identifier.')
  parser.add_argument('-T','--timestamp',
    metavar='TIME',
    dest='timestamp',
    default=TIMESTAMP,
    help='Value to use as the timestamp. Default is the current time. Pretty much any string is OK. Specify a timestamp value or a ##pragma name.')
  parser.add_argument('-c','--chromosomes',
    metavar='CHROMS',
    dest='chromosomes',
    help='Chromosomes. Specifies the chromosome names, preferred order, and lengths. By default, preferred chromosome order is the order they appear in the GFF3 data, and chromosome lengths are set to the max coordinate of their features. This parameter allows you to override both. CHROMS is either a pragma name (typically "##sequence-region" or a string of the form "chr:length,chr:length:...". In the former case, each matching pragma row is parsed. The format of each pragma value is expected to be "chr start end" . ')
  parser.add_argument('-d','--outputDirectory',
    metavar='DIR',
    dest='outputDir',
    default='./output',
    help='Output directory.')
  parser.add_argument('-k','--chunkSize',
    type=int,
    metavar='K',
    dest='chunkSize',
    default=1,
    help='Transcript file chunk size. 0 = everything in one file (no chunking); 1 = one chunk per chromosome; >1 = chunk by chromosome and by start position / chunkSize.')
  parser.add_argument('--sample',
    action='store_true',
    dest='sample',
    default=False,
    help='Sample output.')
  return parser.parse_args()
#
if __name__ == '__main__':
  args = getArgs()
  Gff3Importer(sys.stdin, args).main()

