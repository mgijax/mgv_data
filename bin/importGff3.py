#
#
# importGff3.py
#
# usage:
#    python importGff3.py < MyGenomeAnnotations.gff3
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
    self.fileNameTemplate = '%04d'
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

    self.pragmas = self.datastream.next()
    self.genomeInfo['name'] = self.getOpt('genome')
    self.genomeInfo['timestamp'] = self.getOpt('timestamp')
    #
    chrs = []
    copt = self.getOpt('chromosomes')
    if copt:
      if type(copt) is types.StringType:
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
      'name': self.getOpt('track'),
      'chunkSize': self.opts.chunkSize
    }]
    self.outputDir = os.path.join(self.opts.outputDir, self.genomeInfo['name'].lower().replace('/',''))
    
  def log(self, s):
    sys.stderr.write(s)

  def ensureDirectory (self, d):
    if not os.path.exists(d):
      os.makedirs(d)

  def getFileName(self, chr, blk):
    return os.path.join(self.outputDir, self.getOpt('track'), chr, self.fileNameTemplate % blk)

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

  def writeGrpToBlk(self, grp, chr, blk):
    fname = self.getFileName(chr, blk)
    fd = self.getFileHandle(fname)
    for f in grp:
      if 'exons' in f[8]:
        fd.write(self.formatFeature(f))
        self.outputFiles[fname] += 1

  def writeGrp (self, grp) :
    startBlk = 0 if not self.opts.chunkSize else min(map(lambda f: f[3], grp)) / self.opts.chunkSize
    endBlk = 0 if not self.opts.chunkSize else max(map(lambda f: f[4], grp)) / self.opts.chunkSize
    chr = grp[0][0]
    # if endBlk > startBlk:
    #   self.log('%d .. %d: %s\n' % (startBlk, endBlk, str(grp[0])))
    for blk in range(startBlk, endBlk+1):
      self.writeGrpToBlk(grp, chr, blk)

  def writeTopLevelFeatures(self, lst):
    for f in lst:
      self.tlFile.write(self.formatFeature(f))

  # Processes a gene model group. For top level features, adds a tCount attribute that contains the
  # number of transcripts for that gene. For transcripts, adds an exons attribute whose value is a list 
  # of exon coordinates. Each exon is encoded as 'offset_length', where offset is equal to
  # exon.start - transcript.start, and length is equal to exon.end - exon.start + 1.
  # Also keeps track of all chromosomes seen in the file and the max coordinate seen on each.
  def processGrp(self, grp):
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
    toplevel = []
    for f in grp:
      fid = f[8]['ID']
      if 'Parent' not in f[8]:
        # top level feature
        f[8]['tCount'] = len(pid2kids.get(fid,[]))
        toplevel.append(f)
      elif fid in pid2kids:
        # mid level feature
        exons = pid2kids[f[8]['ID']]
        eexons = f[8].setdefault('exons',[])
        for e in exons:
          offset = e[3] - f[3]
          length = e[4] - f[3] + 1
          eexons.append('%d_%d' % (offset,length))
      else:
        # leaf level feature
        pass
    return toplevel

  def writeGenomeInfo(self):
    fd = open(os.path.join(self.outputDir, 'index.json'), 'w')
    fd.write(json.dumps(self.genomeInfo))
    fd.close()

  def main (self) :
    self.ensureDirectory(self.outputDir)
    self.tlFileName = os.path.join(self.outputDir, 'features.gff3')
    self.tlFile = open(self.tlFileName, 'w')
    for i,grp in enumerate(self.datastream):
      toplevel = self.processGrp(grp)
      self.writeTopLevelFeatures(toplevel)
      self.writeGrp(grp)
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
  parser.add_argument('-t','--trackName',
    metavar='NAME',
    default='models',
    dest='track',
    help='Name of the track.')
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
    default=0,
    help='Chunk size. Set to 0 for no chunking (the default).')
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

