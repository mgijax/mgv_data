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

#
class Gff3Importer:
  def __init__(self, gffSource, genome, track, outputDir, chunkSize, sample):
    self.genome = genome
    self.track = track
    self.gffSource = gffSource
    self.parser = Gff3Parser(self.gffSource, returnHeader=False, returnGroups=True, convertDots=None)
    self.outputDir = outputDir
    self.outputFiles = {}
    self.currFileName = None
    self.currFile = None
    self.chunkSize = chunkSize
    self.sample = sample
    self.fileNameTemplate = 'b%04d.gff3'

  def log(self, s):
    sys.stderr.write(s)

  def ensureDirectory (self, d):
    if not os.path.exists(d):
      os.makedirs(d)

  def getFileName(self, chr, blk):
    gn = self.genome.lower().replace('/','')
    return os.path.join(self.outputDir, gn, chr, self.track, self.fileNameTemplate % blk)

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
      self.currFile.write('[')
      return self.currFile

  def formatFeature(self, f):
    return json.dumps(f)

  def writeGrpToBlk(self, grp, chr, blk):
    fname = self.getFileName(chr, blk)
    fd = self.getFileHandle(fname)
    for f in grp:
      if self.outputFiles[fname]: fd.write(',')
      fd.write(self.formatFeature(f))
      fd.write('\n')
      self.outputFiles[fname] += 1

  def writeGrp (self, grp) :
    startBlk = 0 if not self.chunkSize else min(map(lambda f: f[3], grp)) / self.chunkSize
    endBlk = 0 if not self.chunkSize else max(map(lambda f: f[4], grp)) / self.chunkSize
    chr = grp[0][0]
    # if endBlk > startBlk:
    #   self.log('%d .. %d: %s\n' % (startBlk, endBlk, str(grp[0])))
    for blk in range(startBlk, endBlk+1):
      self.writeGrpToBlk(grp, chr, blk)

  def finalize (self) :
    for fname in self.outputFiles.keys():
      fd = self.getFileHandle(fname)
      fd.write(']')

  def main (self) :
    self.ensureDirectory(self.outputDir)
    for i,grp in enumerate(self.parser.iterate()):
      self.writeGrp(grp)
      if self.sample and i > 1000:
         break
    self.finalize()

#
def getArgs () :
  parser = argparse.ArgumentParser(description='Import genome annotations from a GFF3 file.')
  parser.add_argument('-g','--genomeName',
    metavar='NAME',
    required=True,
    dest='genome',
    help='Name of the genome.')
  parser.add_argument('-t','--trackName',
    metavar='NAME',
    required=True,
    dest='track',
    help='Name of the track.')
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
  Gff3Importer(sys.stdin, genome=args.genome, track=args.track, chunkSize=args.chunkSize, outputDir=args.outputDir, sample=args.sample).main()

